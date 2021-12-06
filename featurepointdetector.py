import math
from featurepoint import FeaturePoint
from django.contrib.gis.geos import fromfile, Point, Polygon

DPI_THRESHOLD = 7
NEIGHBOR_TOTAL_VOTES = 4


def get_id(_id, total):
    if _id > 0:
        _id = _id % total
    else:
        _id = (total + _id) % total
    return _id


def get_vertex_id(point, vertexes):
    for i in range(len(vertexes)):
        if vertexes[i][0] == point[0] and vertexes[i][1] == point[1]:
            return i
    return -1


def get_dpi(polygon, point, distance):
    vertexes = polygon[0][:-1]
    idx = vertexes.index((point[0], point[1]))

    total = 0
    for i in range(math.ceil(len(vertexes) / 2)):
        _n = Point(vertexes[(idx + i) % len(vertexes)]).distance(point)
        _p = Point(vertexes[(idx - i) % len(vertexes)]).distance(point)
        if _n > distance and _p > distance:
            return total

        if _n <= distance:
            total += 1
        if _p <= distance:
            total += 1
    return total


def get_previous(points, actual, positions):
    pos = (len(points) + (actual - positions)) % len(points)
    return Point(points[pos][0], points[pos][1])


def get_next(points, actual, positions):
    pos = (actual + positions) % len(points)
    return Point(points[pos][0], points[pos][1])


class FeaturePointDetector:
    """
    Class based on the implementation of FeaturePointDetector.java by Luı́s Carlos in his thesis "Morphing techniques in
    spatio-temporal database" ("Aplicação de técnicas de morphing em bases de dados espacio-temporais")
    """

    def __init__(self, max_angle, min_angle, d_min, d_max):
        self.max_angle = max_angle
        self.min_angle = min_angle
        self.d_min = d_min
        self.d_max = d_max

        self.dpp = []

    def compute_dpp(self, polygon, distance):
        for p in polygon[0]:
            p = Point(p)
            self.dpp.append(get_dpi(polygon, p, distance))

    def neighbor_vote(self, idx):
        _dpp_sum = 0
        # _dpp_sum = self.dpp[idx]
        for i in range(math.ceil(NEIGHBOR_TOTAL_VOTES / 2)):
            _dpp_sum += self.dpp[(idx + 1) % len(self.dpp)]
            _dpp_sum += self.dpp[(idx - 1) % len(self.dpp)]
        _dpp = (_dpp_sum / NEIGHBOR_TOTAL_VOTES) * .4
        return _dpp + .6 * self.dpp[idx]    # _dpp_sum / (NEIGHBOR_TOTAL_VOTES + 1)

    def get_feature_points(self, polygon):
        vertexes = polygon[0][:-1]  # disregard the last point equal to the first due to the nature of the wkt
        candidates = []
        if len(vertexes) < 2:
            return candidates
        avg_dist = polygon.length / len(polygon[0])
        self.compute_dpp(polygon, 5 * avg_dist)

        # create a list of candidates to feature points
        for i in range(len(vertexes)):
            point = Point(vertexes[i][0], vertexes[i][1])
            feature_point = FeaturePoint(point, 180)
            positions = 1
            p_dist, n_dist = 0, 0

            idx = vertexes.index((point[0], point[1]))
            _dpp = self.neighbor_vote(idx)

            while positions < len(vertexes):
                # get the previous and next vertex of polygon
                _prev = get_previous(vertexes, i, positions)
                _next = get_next(vertexes, i, positions)
                positions += 1

                # verify the d_min condition
                p_dist += point.distance(_prev)
                n_dist += point.distance(_next)
                if self.d_min > (p_dist + n_dist) / 2:
                    continue

                # Verify the d_max condition (one vertex is always tested)
                if positions > 2 and self.d_max > (p_dist + n_dist) / 2:
                    break

                # calculate the angle and check if it's in range
                angle = feature_point.calculate_angle(_prev, _next)
                if self.min_angle < angle < self.max_angle:
                    if angle < feature_point.angle:
                        feature_point.set_angle(angle, _prev, _next)

            # If the vertex is a candidate adds it to the list
            if feature_point.angle < self.max_angle:
                candidates.append(feature_point)

        # print(f"Found candidates: {len(candidates)}/{len(vertexes)} points")

        # filter candidates
        feature_points = []
        for i in range(len(candidates)):
            fp = candidates[i]
            feature_points.append(fp)
            positions = 1
            cont, p_dist, n_dist = True, 0, 0

            while cont and positions < len(candidates):
                # obtain the previous and next candidate
                _prev = candidates[(len(candidates) + (i - positions)) % len(candidates)]
                _next = candidates[(i + positions) % len(candidates)]
                positions += 1
                cont = False

                # compare candidates with the others
                p_dist += fp.point.distance(_prev.point)
                if p_dist < self.d_max:
                    if fp.angle > _prev.angle:
                        feature_points.remove(fp)
                        break
                    cont = True

                n_dist += fp.point.distance(_next.point)
                if n_dist < self.d_max:
                    if fp.angle > _next.angle:
                        feature_points.remove(fp)
                        break
                    cont = True

        return feature_points

    def get_region_of_support_points(self, polygon, fp):
        f_points = self.get_feature_points(polygon)
        vertexes = polygon[0][:-1]

        point_fp_id = get_vertex_id(fp.point, vertexes)
        if point_fp_id < -1:
            print("Error in getRegionOfSupportPoints: invalid Feature Point")

        if len(f_points) < 3:
            ros = []
            if len(f_points) == 0:
                print("Error in getRegionOfSupportPoints: nr of feature points is 0")
            elif len(f_points) == 1:
                points_start = point_fp_id - (len(vertexes) / 2)
                points_end = point_fp_id + (len(vertexes) / 2)

                for i in range(math.floor(points_start), math.floor(points_end)):
                    _id = get_id(i, len(vertexes))
                    if _id == point_fp_id:
                        ros.append(fp)
                    else:
                        ros.append(vertexes[_id])
            elif len(f_points) == 2:
                fp_id = f_points.index(fp)
                if fp_id == 0:
                    # fp_id < prev_id
                    _prev = f_points[1]

                    # get rol
                    points_start = get_vertex_id(_prev.point, vertexes)
                    dist = points_start - point_fp_id
                    if dist > len(vertexes) / 2:
                        dist = len(vertexes) / 2
                    for i in range(dist, 0, -1):
                        _id = get_id(point_fp_id + i, len(vertexes))
                        ros.append(vertexes[_id])

                    ros.append(fp)  # maybe ros.append(fp.point)

                    # get ror
                    dist = (point_fp_id + len(vertexes) - points_start) % len(vertexes)
                    if dist > len(vertexes) / 2:
                        dist = len(vertexes) / 2
                        for i in range(1, math.floor(dist + 1)):
                            _id = get_id(point_fp_id - i, len(vertexes))
                            ros.append(vertexes[_id])

                else:
                    # fp_id > prev_id
                    _prev = f_points[0]

                    # get rol
                    points_start = get_vertex_id(_prev.point, vertexes)
                    dist = point_fp_id - points_start
                    if dist > len(vertexes) / 2:
                        dist = len(vertexes) / 2
                    for i in range(dist, 0, -1):
                        _id = get_id(point_fp_id - i, len(vertexes))
                        ros.append(vertexes[_id])

                    ros.append(fp)  # maybe ros.append(fp.point)

                    # get ror
                    dist = (points_start + len(vertexes) - point_fp_id) % len(vertexes)
                    if dist > len(vertexes) / 2:
                        dist = len(vertexes) / 2
                    for i in range(1, math.floor(dist + 1)):
                        _id = get_id(point_fp_id + 1, len(vertexes))
                        ros.append(vertexes[_id])
            return ros

        max_points = len(vertexes) // len(f_points)
        ros = []
        fp_id = get_vertex_id(fp.point, [p.point for p in f_points])

        # get previous and next fp
        _prev = f_points[get_id(fp_id - 1, len(f_points))]
        _prev_fp_id = get_vertex_id(_prev.point, vertexes)
        _next = f_points[get_id(fp_id + 1, len(f_points))]
        _next_fp_id = get_vertex_id(_next.point, vertexes)

        # limit the region size
        if point_fp_id < _prev_fp_id:
            dist = point_fp_id + len(vertexes) - _prev_fp_id
        else:
            dist = point_fp_id - _prev_fp_id
        if dist > max_points:
            _prev_fp_id = get_id(point_fp_id - max_points, len(vertexes))

        if point_fp_id > _next_fp_id:
            dist = _next_fp_id + len(vertexes) - point_fp_id
        else:
            dist = _next_fp_id - point_fp_id
        if dist > max_points:
            _next_fp_id = get_id(point_fp_id + max_points, len(vertexes))

        # get the ros
        i = 0
        while True:
            _id = get_id(_prev_fp_id + i, len(vertexes))
            i += 1
            if _id == point_fp_id:
                ros.append(fp)
            else:
                ros.append(vertexes[_id])
            if _id == _next_fp_id:
                break
        return ros

    def calculate_feature_points(self, polygon):
        # create a new list of feature points
        ft_points = self.get_feature_points(polygon)
        perimeter = polygon.length

        # calculate and update each feature point proprieties
        for fp in ft_points:
            fp_neighbours = self.get_region_of_support_points(polygon, fp)  # ROS of fp
            fp.calculate_characteristics(fp_neighbours, perimeter)
        return ft_points
