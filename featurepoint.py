import math
from regionofsupport import RegionOfSupport
from django.contrib.gis.geos import Point


class FeaturePoint:
    """
    Class based on the implementation of FeaturePoint2D.java by Luı́s Carlos in his thesis "Morphing techniques in
    spatio-temporal database" ("Aplicação de técnicas de morphing em bases de dados espacio-temporais")
    """

    def __init__(self, point, angle=0):
        self.point = point                  # geos point
        self.angle = angle
        self.prev, self.next = None, None   # geos points
        self.convex = None
        self.ft_size, self.ft_size_l, self.ft_size_r = 0, 0, 0
        self.ft_variation, self.ft_variation_l, self.ft_variation_r = 0, 0, 0
        self.side_ft_variation = 0
        self.n_val, self.t_val = 0, 0
        self.n_vector, self.t_vector = None, None

    def calculate_angle(self, p_minus, p_plus):
        a = self.point.distance(p_plus)
        b = self.point.distance(p_minus)
        c = p_minus.distance(p_plus)
        if a == 0 or b == 0 or c == 0:
            return 0

        # calculate the angle on (ab) and convert it to degrees
        angle = ((a * a) + (b * b) - (c * c)) / (2 * a * b)
        if angle <= -1:
            angle = -1
        elif angle >= 1:
            angle = 1
        return math.degrees(math.acos(angle))

    def set_angle(self, angle, _prev, _next):
        self.angle = angle
        self.prev = _prev
        self.next = _next

        # check convexness
        bx, by = self.point[0] - _prev[0], self.point[1] - _prev[1]
        cx, cy = _next[0] - _prev[0], _next[1] - _prev[1]
        convexness = bx*cy - by*cx
        if convexness >= 0:
            self.convex = True
        else:
            self.convex = False

    def get_least_opening_triangle_points(self):
        return [self.prev, self, self.next]

    def calculate_characteristics(self, fp_neighbours, perimeter):
        # create the ros
        ros = RegionOfSupport(fp_neighbours, self)

        # set fp properties
        self.n_val, self.n_vector = ros.ev_N, ros.e_N
        self.t_val, self.t_vector = ros.ev_T, ros.e_T

        # calculate the feature variation
        self.ft_variation = self.n_val / (self.n_val + self.t_val)
        if not self.convex:
            self.ft_variation *= -1

        # calculate the feature side variation
        fp_id = fp_neighbours.index(self)
        rol = RegionOfSupport(fp_neighbours[0: fp_id+2], self)
        ror = RegionOfSupport(fp_neighbours[fp_id-1: len(fp_neighbours)], self)
        self.ft_variation_l = rol.ev_N / (rol.ev_N + rol.ev_T)
        self.ft_variation_r = ror.ev_N / (ror.ev_N + ror.ev_T)
        if not self.convex:
            self.ft_variation_l *= -1
            self.ft_variation_r *= -1
        self.side_ft_variation = (self.ft_variation_l + self.ft_variation_r) / 2

        # calculate the side variation
        self.ft_size_r, self.ft_size_l = 0, 0
        _prev = fp_neighbours[0]
        for i in range(1, len(fp_neighbours)):
            if i == fp_id:
                _prev = self
                continue
            p = fp_neighbours[i]
            if isinstance(p, tuple):
                p = Point(p)
            else:
                p = p.point

            if not isinstance(_prev, tuple):
                _prev = _prev.point
            if i < fp_id:
                self.ft_size_l += p.distance(Point(_prev[0], _prev[1]))
            else:
                self.ft_size_r += p.distance(Point(_prev[0], _prev[1]))
            _prev = (p[0], p[1])

        self.ft_size_l /= perimeter
        self.ft_size_r /= perimeter
        self.ft_size = (self.ft_size_l + self.ft_size_r) / 2

    def calculate_similarity(self, fp, ft_variation, ft_side, ft_size):
        if self.ft_size > fp.ft_size:
            _max = self.ft_size
        else:
            _max = fp.ft_size
        cost = ft_variation * abs(self.ft_variation - fp.ft_variation)
        cost += ft_side * (abs(self.ft_variation_l-fp.ft_variation_r)+(abs(self.ft_variation_l-fp.ft_variation_l))) / 2
        cost += ft_size * (abs(self.ft_size_l - fp.ft_size_l)+abs(self.ft_size_r - fp.ft_size_r)) / 2
        return _max * cost

    def calculate_discard_cost(self, ft_variation, ft_side, ft_size):
        cost = ft_variation*abs(self.ft_variation) + ft_side*abs(self.side_ft_variation) + ft_size*abs(self.ft_size)
        return cost * self.ft_size

    def __str__(self):
        return f"({self.point.x},{self.point.y}) "
