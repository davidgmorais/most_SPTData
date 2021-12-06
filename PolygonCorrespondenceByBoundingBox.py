import math

from correspondence import Correspondence
from django.contrib.gis.geos import Point, LineString

from featurepoint import FeaturePoint


class PolygonCorrespondenceByBoundingBox:
    """
    Class used to perform polygon correspondence based on the polygon's bounding box, used as an alternative to the
    default method, which works better for bodies whose area/size (and consequently its bounding box area/size) remain
    relatively the same throughout the whole process.
    """

    def __init__(self, s_bbox, t_bbox, s, t):
        """
        Initialization method ot the class

        :param s_bbox: Bounding box of the source polygon
        :param t_bbox: Bounding box of the target polygon
        :param s: Source polygon
        :param t: Target polygon
        """
        self.s_bbox = s_bbox
        self.t_bbox = t_bbox
        self.s = s
        self.t = t
        self.s_fps = []
        self.t_fps = []

    def get_best_approximation(self, origin_point, target_point):
        """
        Auxiliary function to compute the index of the points whose relative distance to the origin is smallest.

        :param origin_point: Origin point
        :param target_point: Target point

        :returns: The index of the target polygon whose relative distance to the origin point is the best.
        """

        _min = math.inf
        _min_idx = None
        for i in range(len(self.t)):
            origin_dist = self.t[i].distance(origin_point)
            target_dist = self.t[i].distance(target_point)
            dist = origin_dist / (origin_dist + target_dist)
            if abs(target_dist - dist) < _min:
                _min_idx = i
        return _min_idx

    def get_min_projections(self, bbox_edges):
        """
        Auxiliary function to compute the indexes and distances of the points of the source polygons whose distance
        to each edge od the bounding box is the smallest.

        :param bbox_edges: Array with the four edges of the bounding box.

        :returns: Array of indexes of points in the source target whose distance is the smallest to a certain edge, and
        an Array of those distances.
        """

        distances, indexes = [], []
        for edge in bbox_edges:

            min_distance, min_index = math.inf, None
            for vertex_index in range(len(self.s)):
                distances_to_edge = edge.project_normalized(self.s[vertex_index])

                if distances_to_edge < min_distance:
                    min_distance = distances_to_edge
                    min_index = vertex_index

            distances.append(min_distance)
            indexes.append(min_index)

        return indexes, distances

    def get_feature_points_correspondence(self):
        """
        Generates the feature point correspondence of each polygon based on their respective bounding box.

        First it breaks the bounding box in LineString objects representing each edge, and then calculates the minimum
        distance of each point of the source polygon to a certain edge. Given this, it goes on to get the point in the
        target polygon whose relative distance is more approximate to the same edge. Finally, it generates an array of
        Correspondences which will contain one correspondence for each edge of the bounding box.

        :returns: Array of Correspondences
        """

        edges_s_bbox = [LineString(
            [Point(self.s_bbox[0][i][0], self.s_bbox[0][i][1]),
             Point(self.s_bbox[0][i + 1][0], self.s_bbox[0][i + 1][1])]
        ) for i in range(len(self.s_bbox[0]) - 1)]

        edges_t_bbox = [LineString(
            [Point(self.t_bbox[0][i][0], self.t_bbox[0][i][1]),
             Point(self.t_bbox[0][i + 1][0], self.t_bbox[0][i + 1][1])]
        ) for i in range(len(self.t_bbox[0]) - 1)]

        s_indexes, s_distances = self.get_min_projections(edges_s_bbox)

        t_indexes = []
        t_distances = []

        for i in range(len(edges_t_bbox)):
            edge = edges_t_bbox[i]
            min_distance, min_index = math.inf, None
            for vertex_index in range(len(self.t)):
                distances_to_edge = edge.project_normalized(self.t[vertex_index])

                diff = abs(distances_to_edge - s_distances[i])
                if diff < min_distance:
                    min_distance = diff
                    min_index = vertex_index

            t_indexes.append(min_index)
            t_distances.append(edge.project_normalized(self.t[min_index]))

        correspondences = []

        for i in range(len(s_indexes)):
            c = Correspondence()
            c.s_i = i
            c.point_s = FeaturePoint(self.s[s_indexes[i]])
            self.s_fps.append(c.point_s)
            c.t_i = i
            c.point_t = FeaturePoint(self.t[t_indexes[i]])
            self.t_fps.append(c.point_t)
            correspondences.append(c)

        return correspondences
