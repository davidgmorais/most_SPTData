import math

from django.contrib.gis.geos import Point, Polygon


def shift_polygon_direction(polygon):
    """
    Auxiliary function to be used to shift the direction of a polygon WKT, meaning if the polygon vertexes are
    clockwise, this whill return them counter-clockwise, and vice-versa.
    """

    vertexes = [Point(x) if isinstance(x, tuple) else Point(x.point) for x in polygon[0][:-1]]
    vertexes.append(vertexes[0])
    return Polygon(vertexes[::-1])


def get_nr_between(begin, end, length):
    """
    Auxiliary function to calculate the number of points between two feature points.

    :param begin: Beginning point to be considered.
    :param end: Ending point to be considered.
    :param length: Length of the polygon.

    :returns: An integer of the number of points between the beginning point and the end point.
    """
    if end > begin:
        return end - begin - 1
    return length - begin + end - 1


def get_weight(begin, end, vertexes):
    """
    Auxiliary function to compute the weigh of each segment between two feature points.

    :param begin: Beginning point to be considered
    :param end: Ending point to be considered
    :param vertexes: Vertexes of the polygon

    :returns: A vector with the weight of each segment from start point to end point
    """

    if end > begin:
        _range = [index for index in range(begin, end)]
    else:
        _range = [index for index in range(begin, len(vertexes))] + [index for index in range(end)]

    weight = [vertexes[i].distance(vertexes[(i + 1) % len(vertexes)]) for i in _range]
    return [x / sum(weight) for x in weight]


def get_remaining_segment(begin, end, index, vertexes):
    """
    Auxiliary function to compute the remaining points in a segment, given an index of said segment.

    :param begin: Beginning point of the segment
    :param end: Ending point of the segment
    :param index: Index of the actual position on the segment
    :param vertexes: Vertexes of the polygon to which the segment belong

    :returns: An array of the remaining vertexes from index to end point
    """

    if end > begin:
        return vertexes[begin + index: end]
    else:
        if begin + index < len(vertexes):
            return vertexes[begin + index:] + vertexes[:end]
        else:
            return vertexes[(begin + index) % len(vertexes): end]


def get_full_segment(begin, end, vertexes):
    """
    Auxiliary function to compute a full segment between two feature points

    :param begin: Beginning poit of the segment
    :param end: Ending point of the segment
    :param vertexes: Vertexes of the polygon to which the segment belongs.

    :return: An array of vertexes belonging to vertexes, from begin point to end point.
    """

    if end > begin:
        return vertexes[begin:end]
    else:
        return vertexes[begin:] + vertexes[:end]


def add_weighted_points(s_segment, t_segment, s_weights, t_weights):
    """
    Add points to the smallest segment based on the weight vectors so that s_segment and t_segment have the same size.

    Add points to the smallest segment by first calculating the relative weights (percentage) of the bigger segment and
    the proceeding to split the smaller segment based on that percentage weigh vector. The result will be two segments,
    the biggest one is unchanged and the smallest one is expanded so that the segment is the same but with more points.

    :param s_segment: Source segment.
    :param t_segment: Target segment.
    :param s_weights: Weight vector of the source segment, its size should be the same of s_segment.
    :param t_weights: Weight vector of the target segment, its size should be the same of t_segment.

    :returns: Two segments of the same size, being the new source segment and the new target segment, respectively.
    """

    if len(t_segment) == len(s_segment):
        return s_segment[:-1], t_segment[:-1]

    if len(t_segment) > len(s_segment):
        biggest_segment, biggest_weights = t_segment, t_weights
        smaller_segment, smaller_weights = s_segment, s_weights
    else:
        biggest_segment, biggest_weights = s_segment, s_weights
        smaller_segment, smaller_weights = t_segment, t_weights

    relative_weights = [w / sum(biggest_weights) for w in biggest_weights]
    # split the smaller segment based on the relative weight
    semi_segment = [smaller_segment[0]]
    for i in range(len(biggest_segment) - 1):
        a, b = smaller_segment[0], smaller_segment[-1]
        d = sum(relative_weights[:i + 1]) * a.distance(b)
        ratio = d / math.sqrt(math.pow(a[0] - b[0], 2) + math.pow(a[1] - b[1], 2))
        _x = a[0] + (b[0] - a[0]) * ratio
        _y = a[1] + (b[1] - a[1]) * ratio
        semi_segment.append(Point(math.floor(_x), math.floor(_y)))

    semi_segment, biggest_segment = semi_segment[:-1], biggest_segment[:-1]
    assert len(semi_segment) == len(biggest_segment)

    if len(t_segment) > len(s_segment):
        return semi_segment, biggest_segment
    return biggest_segment, semi_segment


class PolygonNormalizer:
    """
    Class representing the polygon normalizer, this is, the class that will make so that two polygons have the
    exact same number of points, while not changing its area or form.
    """

    def __init__(self, source, target, correspondences):
        """
        Initializing method Used to normalize two polygons based on their correspondence array.

        Stats by arranging the points of both polygons (and their correspondences) in clockwise direction, if
        they are not already. The normalizes both polygons so they both have the same number of points, and prepares
        them (using the correspondence array) so that the Nth vertex of the source polygon corresponds to the Nth
        vertex of the target polygon.

        :param source: Source polygon to be normalized
        :param target: Target polygon to be normalized
        :param correspondences: Array of Correspondences tho help in the process
        """

        self.s_polygon = source
        self.t_polygon = target
        self.correspondences = correspondences

        c_line_ring = [c.point_s.point for c in correspondences] + [correspondences[0].point_s.point]
        c = Polygon(c_line_ring)
        if c[0].is_counterclockwise:
            # correspondence shift
            self.correspondences.reverse()
        if self.s_polygon[0].is_counterclockwise:
            # source polygon shift
            self.s_polygon = shift_polygon_direction(self.s_polygon)
        if self.t_polygon[0].is_counterclockwise:
            # target polygon shifted
            self.t_polygon = shift_polygon_direction(self.t_polygon)

        self.s_norm = None
        self.t_norm = None

        self.__normalize_polygon()
        self.__prepare_correspondences()

    def __normalize_polygon(self):
        """
        Private auxiliary function to normalize the polygons.

        For each segment between two feature points, calculates the number of points in between for each source and
        target polygon. The according to the number of points in each one uses different strategies to normalize that
        segment so that they both have the same number of points, building the tow final polygon iteratively, one
        segment at a time. In the end generates two Polygon object based on the computed vertexes.
        """

        s_vertexes = [Point(x) if isinstance(x, tuple) else Point(x.point) for x in self.s_polygon[0][:-1]]
        t_vertexes = [Point(x) if isinstance(x, tuple) else Point(x.point) for x in self.t_polygon[0][:-1]]
        _s_polygon, _t_polygon = [], []

        for c in range(len(self.correspondences)):
            _c1 = self.correspondences[c]                                       # first correspondence
            _c2 = self.correspondences[(c + 1) % len(self.correspondences)]     # second correspondence
            s_segment, t_segment = [], []

            # check the index of the fps under analysis in the polygons and see how many points are between each
            s_fp = [s_vertexes.index(_c1.point_s.point), s_vertexes.index(_c2.point_s.point)]
            t_fp = [t_vertexes.index(_c1.point_t.point), t_vertexes.index(_c2.point_t.point)]

            # nr of points in between the source an target polygon between the selected correspondences
            s_between = get_nr_between(s_fp[0], s_fp[1], len(s_vertexes))
            t_between = get_nr_between(t_fp[0], t_fp[1], len(t_vertexes))

            if s_between == t_between:
                # move points to both polygons without adding new ones
                s_segment = get_full_segment(s_fp[0], s_fp[1], s_vertexes)
                t_segment = get_full_segment(t_fp[0], t_fp[1], t_vertexes)

                _s_polygon += s_segment
                _t_polygon += t_segment
                assert len(_s_polygon) == len(_t_polygon)
                continue

            # calculate the weight of each segment between the two feature points under analysis for both polygons
            s_weight = get_weight(s_fp[0], s_fp[1], s_vertexes)
            t_weight = get_weight(t_fp[0], t_fp[1], t_vertexes)

            if s_between == 0:
                # add all the extra points to the source polygon based on weight vector alone
                s_semi_segment = [s_vertexes[s_vertexes.index(_c1.point_s.point)]]
                s_next = [s_vertexes[(s_vertexes.index(s_semi_segment[0]) + 1) % len(s_vertexes)]]

                t_semi_segment = get_full_segment(t_fp[0], t_fp[1], t_vertexes)
                t_next = [t_vertexes[(t_vertexes.index(t_semi_segment[-1]) + 1) % len(t_vertexes)]]

                s_segment, t_segment = add_weighted_points(
                    s_semi_segment + s_next,
                    t_semi_segment + t_next,
                    s_weight, t_weight
                )

            elif t_between == 0:
                # add all the extra points to the target polygon based on weight vector alone
                s_semi_segment = get_full_segment(s_fp[0], s_fp[1], s_vertexes)
                s_next = [s_vertexes[(s_vertexes.index(s_semi_segment[-1]) + 1) % len(s_vertexes)]]

                t_semi_segment = [t_vertexes[t_vertexes.index(_c1.point_t.point)]]
                t_next = [t_vertexes[(t_vertexes.index(t_semi_segment[0]) + 1) % len(t_vertexes)]]

                s_segment, t_segment = add_weighted_points(
                    s_semi_segment + s_next,
                    t_semi_segment + t_next,
                    s_weight, t_weight
                )

            elif s_between < t_between:
                # add points to the source segment based on the relative weight vector
                s_index, t_index = 0, 0

                while t_index < len(t_weight):
                    s_semi_segment, t_semi_segment = [], []
                    s_relative_weight, t_relative_weight = [], []
                    target_sum = 0

                    while True:
                        if t_index >= len(t_weight):
                            target_sum = math.inf
                        else:
                            target_sum += t_weight[t_index]

                        if len(t_semi_segment) > 0 and target_sum >= s_weight[s_index]:
                            s_semi_segment.append(s_vertexes[(s_fp[0] + s_index) % len(s_vertexes)])
                            s_index += 1

                            # if no more points in source, add all to the target semi_segment
                            if s_index == len(s_weight):
                                t_semi_segment += get_remaining_segment(t_fp[0], t_fp[1], t_index, t_vertexes)
                                t_relative_weight += t_relative_weight[t_index:]
                                t_index = len(t_weight)

                            break

                        t_semi_segment.append(t_vertexes[(t_fp[0] + t_index) % len(t_vertexes)])
                        t_relative_weight.append(t_weight[t_index])
                        t_index += 1

                    s_next = [s_vertexes[(s_vertexes.index(s_semi_segment[-1]) + 1) % len(s_vertexes)]]
                    t_next = [t_vertexes[(t_vertexes.index(t_semi_segment[-1]) + 1) % len(t_vertexes)]]
                    s_semi_segment, t_semi_segment = add_weighted_points(
                        s_semi_segment + s_next,
                        t_semi_segment + t_next,
                        s_relative_weight,
                        t_relative_weight
                    )

                    s_segment += s_semi_segment
                    t_segment += t_semi_segment

            else:
                # add points to the target segment based on the relative weight vector
                s_index, t_index = 0, 0

                while s_index < len(s_weight):
                    s_semi_segment, t_semi_segment = [], []
                    s_relative_weight, t_relative_weight = [], []
                    source_sum = 0

                    while True:
                        if s_index >= len(s_weight):
                            source_sum = math.inf
                        else:
                            source_sum += s_weight[s_index]

                        if len(s_semi_segment) > 0 and source_sum >= t_weight[t_index]:
                            t_semi_segment.append(t_vertexes[(t_fp[0] + t_index) % len(t_vertexes)])
                            t_relative_weight.append(t_weight[t_index])
                            t_index += 1

                            # if no more points in target, add all to the source semi_segment
                            if t_index == len(t_weight):
                                s_semi_segment += get_remaining_segment(s_fp[0], s_fp[1], s_index, s_vertexes)
                                s_relative_weight += s_relative_weight[s_index:]
                                s_index = len(s_weight)

                            break

                        s_semi_segment.append(s_vertexes[(s_fp[0] + s_index) % len(s_vertexes)])
                        s_relative_weight.append(s_weight[s_index])
                        s_index += 1

                    s_next = [s_vertexes[(s_vertexes.index(s_semi_segment[-1]) + 1) % len(s_vertexes)]]
                    t_next = [t_vertexes[(t_vertexes.index(t_semi_segment[-1]) + 1) % len(t_vertexes)]]
                    s_semi_segment, t_semi_segment = add_weighted_points(
                        s_semi_segment + s_next,
                        t_semi_segment + t_next,
                        s_relative_weight,
                        t_relative_weight
                    )

                    s_segment += s_semi_segment
                    t_segment += t_semi_segment

            assert len(s_segment) == len(t_segment)
            _s_polygon += s_segment
            _t_polygon += t_segment

        _s_polygon.append(_s_polygon[0])
        _t_polygon.append(_t_polygon[0])

        assert len(_s_polygon) == len(_t_polygon)
        self.s_norm = Polygon([(p[0], p[1]) for p in _s_polygon])
        self.t_norm = Polygon([(p[0], p[1]) for p in _t_polygon])

    def __prepare_correspondences(self):
        """
        Private auxiliary function to help with rotating the polygons so that their vertexes correspond.

        This is done by using the first correspondence as a starting point and then starting to read points from both
        polygons. We end up with two new polygons with the same points, but to each applied a certain offset
        mandated by the first correspondence.
        """

        s_vertexes = [Point(x) if isinstance(x, tuple) else Point(x.point) for x in self.s_norm[0][:-1]]
        t_vertexes = [Point(x) if isinstance(x, tuple) else Point(x.point) for x in self.t_norm[0][:-1]]
        _c = self.correspondences[0]
        s_shifted, t_shifted = [], []

        s_index = s_vertexes.index(_c.point_s.point)
        t_index = t_vertexes.index(_c.point_t.point)
        assert len(s_vertexes) == len(t_vertexes)

        for i in range(len(t_vertexes)):
            s_shifted.append(s_vertexes[(i + s_index) % len(s_vertexes)])
            t_shifted.append(t_vertexes[(i + t_index) % len(t_vertexes)])

        s_shifted.append(s_shifted[0])
        t_shifted.append(t_shifted[0])

        self.s_norm = Polygon([(p[0], p[1]) for p in s_shifted])
        self.t_norm = Polygon([(p[0], p[1]) for p in t_shifted])

