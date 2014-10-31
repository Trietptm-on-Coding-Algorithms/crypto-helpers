import math

def STR(a):
    # Yes, this is a little bit silly :-)
    for i in range(0, len(a) - 1, 2):
            print((chr(int(a[i:i+2]))), end=' ')

def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def modinv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m

class Coord(object):
    def __init__(self, x, field):
        self.x = x
        self.n = n

    def __add__(self, other):
        if isinstance(other, Coord):
            if self.n != other.n:
                raise Exception('different fields')
            other = other.x

        return Coord((self.x + other) % self.n, self.n)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        if isinstance(other, Coord):
            if self.n != other.n:
                raise Exception('different fields')
            other = other.x

        return Coord((self.x - other) % self.n, self.n)

    def __mul__(self, other):
        if isinstance(other, Coord):
            if self.n != other.n:
                raise Exception('different fields')
            other = other.x

        res = self.x * other
        return Coord(res % self.n, self.n)

    def __rmul__(self, other):
        return self*other

    def __pow__(self, k):
        return Coord(pow(self.x, k, self.n), self.n)

    def __truediv__(self, other):
        if isinstance(other, Coord):
            other = other.x

        return Coord(self.x*modinv(other, self.n), self.n)

    def __eq__(self, other):
        if isinstance(other, Coord):
            if self.n != other.n:
                raise Exception('different fields')
            other = other.x

        return (self.x == other)

    def __ne__(self, other):
        if isinstance(other, Coord):
            if self.n != other.n:
                raise Exception('different fields')
            other = other.x

        return (self.x != other)

    def __str__(self):
        return "%s" % self.x

class Point(object):
    def __init__(self, x, y, n):
        if not isinstance(x, Coord):
            x = Coord(x, n)

        if not isinstance(y, Coord):
            x = Coord(y, n)

        self.x = x
        self.y = y
        self.n = n

    def __str__(self):
        return "(%s, %s)" % (self.x, self.y)

class PointEq(Point):
    def __init__(self, x, y, eq):
        super().__init__(x, y, eq.n)
        self.eq = eq

    def __mul__(self, k):
        return self.eq.mul(self, k)

    def __rmul__(self, other):
        return self*other

class Equation(object):
    def __init__(self, params, field):
        self.a = params[0]
        self.b = params[1]
        self.n = field

    def __getitem__(self, key):
        x = Coord(key[0], self.n)
        y = Coord(key[1], self.n)

        if (y**2 != x**3 + self.a*x + self.b):
            raise Exception('Point is not on curve (%s, %s)' % (y**2, x**3 + x*self.a + self.b))

        return PointEq(x, y, self)

    def add(self, p, q):
        if (p.x == None) and (p.y == None):
            return q

        if (q.x == None) and (q.y == None):
            return self

        if (p.x==q.x) and (p.y==q.y):
            return self.double(p)

        s = (p.y - q.y)/(p.x - q.x)
        rx = (s*s - p.x - q.x)
        ry = (s*(p.x - rx) - p.y)

        return PointEq(rx, ry, eq)

    def double(self, p):
        s = (3*p.x*p.x + self.a)/(2*p.y)
        rx = (s*s - 2*p.x)
        ry = (s*(p.x - rx) - p.y)

        return PointEq(rx, ry, self)
    
    def mul(self, p, k):
        N = p
        R = Point(None, None, self.n)

        for bit in range(0, k.bit_length()+1):
            if (k & (1<<bit)):
                R = self.add(R, N)
            N = self.double(N)

        return R

    def mod(self, n):
        return PointEq(self.x % n, self.y % n, self)

if __name__ == "__main__":
    e = 141597355687225811174313
    d = 87441340171043308346177
    n = 928669833265826932708591
    C = (236857987845294655469221, 12418605208975891779391)

    eq = Equation([0, pow(C[1],2,n) - pow(C[0],3, n)], n)
    c = eq[C[0], C[1]]

    #print eq.double(c) # (196269932798073152664662 : 66403428568175199097197 : 1)
    #print eq.mul(c, 3) # (154214382988446711829623 : 588038293798980333372500 : 1)

    M = c * d
    C2 = M * e

    print(C2, M)
    print(M.x)
    print(M.y)
    STR(str(M.x)), STR(str(M.y))
