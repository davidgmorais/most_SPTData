import concurrent.futures
import getopt
import json
import sys
import re
import math
import os
import time

from django.contrib.gis.geos import Polygon, Point

from featurepointdetector import FeaturePointDetector
from PolygonCorrespondenceByBoundingBox import PolygonCorrespondenceByBoundingBox
from polygoncorrespondence import PolygonCorrespondence
from polygonNormalizer import PolygonNormalizer


DEFAULT_CONFIGS = {
    "max_angle": 160,
    "min_angle": 20,
    "d_min": 2.8,
    "d_max": 1.5,
    "simplification_folder": "/simplification"
}


class Tracker:
    """
    Class used to help with the man_track file information
    """

    def __init__(self, man_track_file):
        """
        Generates an array of entries based on the man_track_file where, accounting to
        https://public.celltrackingchallenge.net/documents/Naming%20and%20file%20content%20conventions.pdf, the first
        column represent the label of an entity, the second and third represent the first and last frame in which they
        appear (respectively), and the last one represents the entity's parent label. Furthermore, it also calculates
        which frame is considered the first so the second columns does not need to necessarily start at 0.

        :param man_track_file: Path to the file man_track.txt representing the acyclic graphic for the video.
        """

        self.entities = []
        f = open(man_track_file, 'r')
        start_frame = math.inf
        for _line in f.readlines():
            label, start, end, parent = _line.strip().split(" ")
            self.entities.append((int(label), int(start), int(end), int(parent)))
            if int(start) < start_frame:
                start_frame = int(start)
        f.close()
        self.starting_frame = start_frame

    def on_frame(self, frame_nr):
        """
        Computes the array of entities present in a given frame.

        :param frame_nr: Number of the frame to be considered.

        :returns: Array of entities in the frame passed as parameter.
        """

        entities_on_frame = []
        for e in self.entities:
            if e[1] <= frame_nr <= e[2]:
                entities_on_frame.append(e)

        return entities_on_frame

    def get_label_by_index_on_frame(self, index, frame):
        """
        Get the label of an entity, based on their index, in any given frame.

        :param index: Index of the entity in the frame
        :param frame: Number of frame to be considered

        :returns: Entity which index and frame match the ones passed as parameter
        """

        entities_on_frame = self.on_frame(self.starting_frame + frame)
        return entities_on_frame[index]

    def get_parent_by_label(self, label):
        """
        Get the label of an entity's parent, based on their label.

        :param label: Label of the entity whose parent entity's label we want

        :returns: The entity's parent label, if it exists, 0 otherwise.
        """
        for entity in self.entities:
            if entity[0] == label:
                return entity[3]
        return 0


def get_configs_from_file(config_filename):
    """
    Loads the configuration dictionary from a file. If a certain key is missing on the configuration file, the default
    value present at DEFAULT_CONFIGURATION for that key is used.

    :param config_filename: Path to the configuration JSON file to be imported

    :returns: Configuration dictionary based on the imported file.
    """

    config_file = open(config_filename, "r")
    configs = json.load(config_file)

    for key in DEFAULT_CONFIGS.keys():
        if key not in configs.keys():
            print(f"[WARNING] {key} configuration was not provided - using default value ({DEFAULT_CONFIGS[key]})")
            configs[key] = DEFAULT_CONFIGS[key]

    return configs


def read_polygon_from_file(filename):
    """
    Read a polygon's WKT from a file and creates a Polygon object using their vertexes. If the WKT is malformed and the
    first point does not match the last, then the first point is copied to the end of the vertexes array (this is why
    this function is used instead of the fromfile() built-in method).

    :param filename: Path to the file containing the WKT to be imported

    :returns: Polygon object based on the vertexes present in the WKT file passed as a parameter.
    """

    f = open(filename, 'r')
    wkt = f.readline().split("((")[-1].split("))")[0]
    coordinates = [coordinate for coordinate in wkt.split(", ")]
    coordinates = [(float(c.split()[0]), float(c.split()[1])) for c in coordinates]
    if coordinates[-1] != coordinates[0]:
        coordinates.append(coordinates[0])

    return Polygon(coordinates)


def rotate(polygon, theta, center):
    """
    Auxiliary function to help on the creation of a polygon's bounding box.

    Rotates a polygon by Rotating all its points by a certain angle relative to a center and appends them to a new
    array, using that array to create a new Polygon object in order to apply the bounding box method to it.

    :param polygon: Polygon to be rotated.
    :param theta: Angle of the rotation.
    :param center: Center the rotation.

    :returns: Polygon object representing the rotated polygon.
    """

    cos_ang, sin_ang = math.cos(theta), math.sin(theta)
    points = polygon[0][:-1]

    new_polygon = []
    for p in points:
        x, y = p[0], p[1]
        tx, ty = x - center[0], y - center[1]
        new_x = (tx * cos_ang + ty * sin_ang) + center[0]
        new_y = (-tx * sin_ang + ty * cos_ang) + center[1]
        new_polygon.append(Point(new_x, new_y))
    new_polygon.append(new_polygon[0])
    return Polygon(new_polygon)


def create_bounding_box(_p):
    """
    Creates a bounding box of a polygon, first rotating it so that the points with bigger and smaller x are in the
    horizontal, and then applies the built-in method envelop of a Polygon to it. Afterwards it rotates the created
    bounding box in the opposite direction to math the angle of the original polygon.

    :param _p: Polygon to be rotated.

    :returns: Bounding box of the polygon passed as a parameter.
    """

    s_vertex = _p[0][:-1]

    s_vertex.sort(key=lambda p: p[0])
    _x = s_vertex[::len(s_vertex) - 1]
    _x_angle = math.atan2(_x[1][1] - _x[0][1], _x[1][0] - _x[0][0])

    _p_rotated = rotate(_p, _x_angle, _p.centroid)
    _p_rotated_bbox = _p_rotated.envelope
    _p_bbox = rotate(_p_rotated_bbox, -_x_angle, _p_rotated.centroid)
    return _p_bbox


def vcp(index, s, t, s_simplification, t_simplification, bbox_method, max_angle, min_angle, d_min, d_max):
    """
    Solves the Vertex Correspondence Problem for one pair of polygons (source and target) using the max and min angles
    and distances presented in the configuration dictionary and (optionally) simplified version of the WKTs for the
    source and target polygon.

    The default method by which this problem is solved is based on the solution proposed and implemented by Luı́s Carlos
    Marques Paulo in his thesis "Morphing techniques in spatio-temporal databases" ("Aplicação de técnicas de morphing
    em bases de dados espacio-temporais"), but an alternative method was created in which the bounding boxes of the
    polygons are used to create a correspondence between vertexes. We observed that this optional method was better
    suited for entities whose area did not suffer great alterations during the video (such as cells and icebergs) but
    still was lackluster for fires when compared to the default method.

    :param index: Index of the frame under analysis.
    :param s: Source Polygon object.
    :param t: Target Polygon object.
    :param s_simplification: Polygon object representing the Simplified version of the source polygon.
    :param t_simplification: Polygon object representing the Simplified version of the target polygon.
    :param bbox_method: Flag to identify the method used (True uses the bounding box method, False uses the default).
    :param max_angle: Max angle used in feature point detection.
    :param min_angle: Min angle used in feature point detection.
    :param d_max: Max distance used in feature point detection.
    :param d_min: Min distance used in feature point detection.

    :returns: The index of the frame under analysis, as well as both the source and target polygons who underwent
    normalization and the correspondences generated by the method.
    """

    if bbox_method:
        _p1_bbox = create_bounding_box(s)
        _p2_bbox = create_bounding_box(s)
        pc = PolygonCorrespondenceByBoundingBox(_p1_bbox, _p2_bbox,
                                                [Point(x) for x in s[0][:-1]],
                                                [Point(x) for x in t[0][:-1]])
        correspondences = pc.get_feature_points_correspondence()

    else:
        # Feature point extraction
        n, m = len(s[0]), len(t[0])
        s_avg_dist = s.length / n
        t_avg_dist = t.length / m

        fp_detector = FeaturePointDetector(
            max_angle, min_angle, s_avg_dist / d_min, s_avg_dist * d_max
        )
        if s_simplification:
            s_feature_points = fp_detector.calculate_feature_points(s_simplification)
        else:
            s_feature_points = fp_detector.calculate_feature_points(s)

        fp_detector = FeaturePointDetector(
            max_angle, min_angle, t_avg_dist / d_min, t_avg_dist * d_max
        )
        if t_simplification:
            t_feature_points = fp_detector.calculate_feature_points(t_simplification)
        else:
            t_feature_points = fp_detector.calculate_feature_points(t)

        print(f"\t- Generated feature points {index}-{index + 1}")

        # compute correspondences
        pc1 = PolygonCorrespondence(s_feature_points, t_feature_points, .3, .3, .3, 1)
        correspondences = pc1.get_feature_point_correspondences(2)
        print(f"\t- Generated correspondence {index}-{index + 1}")

    normalizer = PolygonNormalizer(s, t, correspondences)
    s_normalized = normalizer.s_norm
    t_normalized = normalizer.t_norm

    return index, s_normalized, t_normalized, correspondences


def vcp_man_track(working_dir, man_track_file, config_filename=None, bbox_method=False, output_filename=None):
    print("Checking configs and folders on working directory...")
    if config_filename:
        configs = get_configs_from_file(config_filename)
    else:
        configs = DEFAULT_CONFIGS

    directories = sorted([int(_d) for _d in os.listdir(working_dir) if os.path.isdir(working_dir + "/" + _d)])
    tracker = Tracker(man_track_file)
    correspondences_dict = dict()

    # for loop
    for i in range(len(directories) - 1):
        source_frame = i
        target_frame = i + 1

        source_dir = working_dir + "/" + str(source_frame)
        target_dir = working_dir + "/" + str(target_frame)

        source_polygons_filename = [label.split(".wkt")[0] for label in os.listdir(source_dir)]
        target_polygons_filename = [label.split(".wkt")[0] for label in os.listdir(target_dir)]

        s_frame = tracker.on_frame(tracker.starting_frame + i)
        t_frame = tracker.on_frame(tracker.starting_frame + i + 1)
        assert len(source_polygons_filename) <= len(s_frame)
        assert len(target_polygons_filename) <= len(t_frame)

        s_frame = s_frame[:len(source_polygons_filename)]
        t_frame = t_frame[:len(target_polygons_filename)]

        print("Reading WKT files...")
        source_polygons = []
        target_polygons = []
        for j in range(len(s_frame)):
            s_polygon = read_polygon_from_file(source_dir + "/" + source_polygons_filename[j] + ".wkt")
            source_polygons.append(s_polygon)

        for j in range(len(t_frame)):
            t_polygon = read_polygon_from_file(target_dir + "/" + target_polygons_filename[j] + ".wkt")
            target_polygons.append(t_polygon)

        print(f"*** Analysing correspondence: {source_frame} - {target_frame} ***")
        print("Generating Feature points and correspondences...")
        start = time.time()
        with concurrent.futures.ProcessPoolExecutor(max_workers=6) as executor:
            _future = {
                executor.submit(vcp, i, source_polygons[i], target_polygons[i + 1], None, None, bbox_method,
                                configs['max_angle'], configs['min_angle'], configs['d_min'], configs['d_max']): i
                for i in range(min(len(source_polygons), len(target_polygons)) - 1)
            }

            done, _ = concurrent.futures.wait(_future, timeout=None, return_when=concurrent.futures.ALL_COMPLETED)
            for future in done:
                index, s_normalized, t_normalized, correspondences = future.result()

                if source_frame in correspondences_dict.keys():
                    correspondences_dict[source_frame][index] = [s_normalized, t_normalized, correspondences]
                else:
                    correspondences_dict[source_frame] = {index: [s_normalized, t_normalized, correspondences]}

        end = time.time()
        print("Process pool finished in {:.3f}s".format(end - start))

    print("Generating GeoJson file...")
    features = []
    for frame in sorted(correspondences_dict.keys()):
        geometries = []

        for polygon in sorted(correspondences_dict[frame].keys()):
            s_label = tracker.get_label_by_index_on_frame(polygon, frame)
            s_geometry = {
                "type": "Polygon",
                "properties": {
                    "polygon-type": "Source",
                    "id": polygon,
                    "label": s_label[0],
                    "parent": tracker.get_parent_by_label(s_label[0])
                },
                "coordinates": [correspondences_dict[frame][polygon][0][0][:]]
            }
            geometries.append(s_geometry)

            t_label = tracker.get_label_by_index_on_frame(polygon, frame)
            t_geometry = {
                "type": "Polygon",
                "properties": {
                    "polygon-type": "Target",
                    "id": polygon,
                    "label": t_label[0],
                    "parent": tracker.get_parent_by_label(t_label[0])
                },
                "coordinates": [correspondences_dict[frame][polygon][1][0][:]]
            }
            geometries.append(t_geometry)

        feature = {
            "type": "Feature",
            "properties": {
                "order": frame,
            },
            "geometry": {
                "type": "GeometryCollection",
                "geometries": geometries
            }
        }
        features.append(feature)

    geo_json = {
        "type": "FeatureCollection",
        "properties": {
            "total frames": len(correspondences_dict.keys()) + 1,
        },
        "features": features
    }
    geo_json = json.dumps(geo_json, indent=4)

    output_filename = working_dir + '/' + output_filename
    jsonFile = open(output_filename, "w")
    jsonFile.write(geo_json)
    jsonFile.close()
    print(f"Output stored in {output_filename}")


def vcp_multiple(working_dir, config_filename=None, bbox_method=False, output_filename=None):
    print("Checking configs and folders on working directory...")
    if config_filename:
        configs = get_configs_from_file(config_filename)
    else:
        configs = DEFAULT_CONFIGS

    simplification_folder = False
    if os.path.exists(working_dir + str(configs['simplification_folder'])):
        simplification_folder = True

    print("Reading WKT files...")
    wkt_list = [wkt for wkt in os.listdir(working_dir)
                if os.path.isfile(working_dir + "/" + wkt) and wkt.split(".")[-1] == "wkt"]
    wkt_list.sort(key=lambda wkt: int(wkt.split(".wkt")[0]))
    correspondences_dict = dict()

    polygons, simplifications = [], []
    for i in range(len(wkt_list)):
        polygons.append(read_polygon_from_file(working_dir + "/" + wkt_list[i]))
        if simplification_folder and os.path.isfile(
                working_dir + str(configs['simplification_folder']) + "/" + wkt_list[i]):
            simplifications.append(
                read_polygon_from_file(working_dir + str(configs['simplification_folder']) + "/" + wkt_list[i])
            )
        else:
            simplifications.append(None)

    print("Generating Feature points and correspondences...")
    start = time.time()
    with concurrent.futures.ProcessPoolExecutor(max_workers=3) as executor:
        _future = {
            executor.submit(vcp, i, polygons[i], polygons[i + 1], simplifications[i], simplifications[i + 1],
                            bbox_method,
                            configs['max_angle'], configs['min_angle'], configs['d_min'], configs['d_max']): i
            for i in range(len(wkt_list) - 1)
        }

        done, _ = concurrent.futures.wait(_future, timeout=None, return_when=concurrent.futures.ALL_COMPLETED)
        for future in done:
            index, s_normalized, t_normalized, correspondences = future.result()
            correspondences_dict[index] = [s_normalized, t_normalized, correspondences]

    end = time.time()
    print("Process pool finished in {:.3f}s".format(end - start))

    if output_filename:
        output_filename = working_dir + '/' + output_filename
        print("Generating GeoJson file...")
        features = []
        for index in sorted(correspondences_dict.keys()):
            geometries = []
            s_geometry = {
                "type": "Polygon",
                "properties": {"polygon-type": "Source"},
                "coordinates": [correspondences_dict[index][0][0][:]],
            }
            geometries.append(s_geometry)
            t_geometry = {
                "type": "Polygon",
                "properties": {"polygon-type": "Target"},
                "coordinates": [correspondences_dict[index][1][0][:]],
            }
            geometries.append(t_geometry)

            feature = {
                "type": "Feature",
                "properties": {
                    "order": index,
                },
                "geometry": {
                    "type": "GeometryCollection",
                    "geometries": geometries
                }
            }
            features.append(feature)

        geo_json = {
            "type": "FeatureCollection",
            "properties": {
                "total frames": len(correspondences_dict.keys()) + 1,
            },
            "features": features
        }
        geo_json = json.dumps(geo_json, indent=4)

        jsonFile = open(output_filename, "w")
        jsonFile.write(geo_json)
        jsonFile.close()
        print(f"Output stored in {output_filename}")

    else:
        print("\n*** Correspondences generated ***")
        for index in sorted(correspondences_dict.keys()):
            print()
            print(f"From {index} to {index + 1}:")
            for c in correspondences_dict[index][2]:
                print(f"\t{c.point_s.point} -> {c.point_t.point}")


def vcp_single(s_file, t_file, s_simplification_file=None, t_simplification_file=None, config_files=None, bbox=False):
    source = read_polygon_from_file(s_file)
    target = read_polygon_from_file(t_file)

    s_simplification, t_simplification = None, None
    if s_simplification_file and os.path.isfile(s_simplification_file):
        s_simplification = read_polygon_from_file(s_simplification_file)
    if t_simplification_file and os.path.isfile(t_simplification_file):
        t_simplification = read_polygon_from_file(t_simplification_file)

    if config_files:
        configs = get_configs_from_file(config_files)
    else:
        configs = DEFAULT_CONFIGS

    index, s_normalized, t_normalized, correspondences = vcp(
        0, source, target, s_simplification, t_simplification, bbox,
        configs['max_angle'], configs['min_angle'], configs['d_min'], configs['d_max']
    )
    print("Generated correspondence:")
    print("Source: ", s_normalized)
    print("Target: ", t_normalized)


if __name__ == "__main__":
    _dir = None
    source_filename, target_filename = None, None
    man_track_filename = None

    if len(sys.argv) <= 1:
        print("[ERROR] Wrong parameters, check help with -h")

    elif len(sys.argv) == 2 or re.match('^-[a-zA-Z]$', sys.argv[2]):
        if sys.argv[1] != '-h':
            _dir = sys.argv[1]
    else:
        if os.path.isfile(sys.argv[1]) and os.path.isfile(sys.argv[2]):
            print("single processing")
            source_filename = sys.argv[1]
            target_filename = sys.argv[2]
        elif os.path.isdir(sys.argv[1]) and os.path.isfile(sys.argv[2]):
            print("man_track processing")
            _dir = sys.argv[1]
            man_track_filename = sys.argv[2]
        else:
            print("[ERROR] Wrong parameters, check help with -h")

    # optional handling
    if _dir:
        if man_track_filename:
            argv = 3
        else:
            argv = 2
        options = "bo:c:h"
        long_options = ["bbox", "output-file", "config-file", "help"]
    elif source_filename and target_filename:
        argv = 3
        options = "s:t:c:bh"
        long_options = ["source-simplification", "target-simplification", "config-file", "bbox", "help"]
    else:
        argv = 1
        options = "h"
        long_options = ["help"]

    opts, args = None, None
    try:
        opts, args = getopt.getopt(sys.argv[argv:], options, long_options)
    except getopt.GetoptError:
        print("[ERROR] Error parsing options")

    _s_simplification, _t_simplification = None, None
    _bbox_method = False
    _output_filename, _config_file = None, None

    for opt, arg in opts:
        if opt in ["-h", "--help"]:
            usage = open("usage.txt", "r")
            for line in usage.readlines():
                print(line, end='')
            usage.close()
            exit(0)
        elif opt in ["-s", "--source-simplification"] and not _dir:
            _s_simplification = arg
        elif opt in ["-t", "--target-simplification"] and not _dir:
            _t_simplification = arg
        elif opt in ["-b", "--bbox"]:
            _bbox_method = True
        elif opt in ["-o", "--output-file"]:
            _output_filename = arg
        elif opt in ["-c", "--config-file"]:
            _config_file = arg
            if _config_file.split(".")[-1] != "json":
                print("[ERROR] Please provide a valid JSON configuration file.")
                exit(1)
        else:
            print("[ERROR] Wrong usage, please check help with -h")

    if _dir:
        if not os.path.isdir(_dir):
            print("[ERROR]", _dir, "is not a directory")
            exit(1)
        if _dir[-1] == "/":
            _dir = _dir[:-1]

        if man_track_filename:
            vcp_man_track(_dir, man_track_filename, _config_file, _bbox_method, _output_filename)
        else:
            vcp_multiple(_dir, _config_file, _bbox_method, _output_filename)

    else:
        vcp_single(source_filename, target_filename, _s_simplification, _t_simplification, _config_file, _bbox_method)
