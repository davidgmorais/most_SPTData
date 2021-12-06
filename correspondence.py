class Correspondence:
    """
    Class based on the implementation of Correspondence.java by Luı́s Carlos in his thesis "Morphing techniques in
    spatio-temporal database" ("Aplicação de técnicas de morphing em bases de dados espacio-temporais")
    """

    def __init__(self):
        self.point_s = None
        self.s_i = 0
        self.point_t = None
        self.t_i = 0
        self.cost = 0

    def __str__(self):
        return f"[{self.s_i}]({self.point_s}) -> [{self.t_i}]({self.point_t});"

