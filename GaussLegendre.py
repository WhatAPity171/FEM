import math
import inspect
class GaussLegendreParams:
    def __init__(self, n):
        if n not in [1, 2, 3, 4]:
            raise ValueError("Obsługiwane są tylko wartości n = 1, 2, 3 lub 4")

        self.n = n
        self.points, self.weights = self._get_points_and_weights(n)

    def _get_points_and_weights(self, n):
        if n == 1:
            points = [0.0]
            weights = [2.0]
        elif n == 2:
            points = [-1/math.sqrt(3), 1/math.sqrt(3)]
            weights = [1.0, 1.0]
        elif n == 3:
            points = [-math.sqrt(3/5), 0.0, math.sqrt(3/5)]
            weights = [5/9, 8/9, 5/9]
        elif n == 4:
            points = [
                -math.sqrt((3 + 2*math.sqrt(6/5)) / 7),
                -math.sqrt((3 - 2*math.sqrt(6/5)) / 7),
                math.sqrt((3 - 2*math.sqrt(6/5)) / 7),
                math.sqrt((3 + 2*math.sqrt(6/5)) / 7)
            ]
            weights = [
                (18 - math.sqrt(30)) / 36,
                (18 + math.sqrt(30)) / 36,
                (18 + math.sqrt(30)) / 36,
                (18 - math.sqrt(30)) / 36
            ]
        return points, weights
    
def f1(x : float) -> float: 
    return 5*x*x + 3*x +6

def f2(x : float, y : float) -> float:
    return 5*x*x*y*y + 3*x*y + 6  

def GaussIntegration(params : GaussLegendreParams, function):
    result = 0.0
    n_params = len(inspect.signature(function).parameters)
    if n_params == 1:
        for xi,wi in zip(params.points, params.weights):
            result += function(xi) * wi
    if n_params == 2:
        for xi,wxi in zip(params.points, params.weights):
            for yi,wyi in zip(params.points, params.weights):
                result += function(xi,yi) * wxi * wyi
    return result


if __name__ == "__main__":
    params = GaussLegendreParams(2)
    print(GaussIntegration(params, f2)) 