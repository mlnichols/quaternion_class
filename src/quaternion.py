from __future__ import (absolute_import,division,print_function,unicode_literals)

from math import exp,log,cos,sin,acos
class quaternion:
    def __init__(self,w=0.,x=0.,y=0.,z=0.):
        self.w = w
        self.x = x
        self.y = y
        self.z = z
    def __repr__(self):
        return "%g%+gi%+gj%+gk" % (self.w, self.x, self.y, self.z)
    def __neg__(self):
        return quaternion(-self.w,-self.x,-self.y,-self.z)
    def __pos__(self):
        return quaternion(self.w,self.x,self.y,self.z)
    def __add__(self,other):
        if isinstance(other,quaternion):
            return quaternion(self.w+other.w,self.x+other.x,self.y+other.y,self.z+other.z)
        else:
            raise TypeError("Quaternions can only be added to other quaternions to avoid ambiguity.")
    def __radd__(self,other):
        raise TypeError("Quaternions can only be added to other quaternions to avoid ambiguity.")
    def __sub__(self,other):
        if isinstance(other,quaternion):
            return self + other._neg__()
        else:
            raise TypeError("Quaternions can only be substracted by other quaternions to avoid amibiguity.")
    def __rsub__(self,other):
        raise TypeError("Quaternions can only be substracted by other quaternions to avoid amibiguity.")
    def __mul__(self,other):
        if isinstance(other,complex):
            raise TypeError("Multiplication by a complex number is disallowed due to ambiguity.")
        elif isinstance(other,quaternion):
            return quaternion(
                self.w*other.w - self.x*other.x - self.y*other.y - self.z*other.z,
                self.w*other.x + self.x*other.w + self.y*other.z - self.z*other.y,
                self.w*other.y - self.x*other.z + self.y*other.w + self.z*other.x,
                self.w*other.z + self.x*other.y - self.y*other.x + self.z*other.w
                )
        else:
            return quaternion(self.w*other,self.x*other,self.y*other,self.z*other)
    def __rmul__(self,other):
        if isinstance(other,complex):
            raise TypeError("Multiplication by a complex number is disallowed due to ambiguity.")
        else:
            return quaternion(self.w*other,self.x*other,self.y*other,self.z*other)
    def __conjugate__(self):
        return quaternion(self.w,-self.x,-self.y,-self.z)
    def __abs__(self):
        return (self.w*self.w + self.x*self.x + self.y * self.y + self.z*self.z)**0.5
    def __div__(self,other): #Redudant due to future import. Kept for readability
        try:
            fother = float(other)
        except:
            raise TypeError("Quaternions can only be divided by a real number.")
        return quaternion(self.w/fother,self.x/fother,self.y/fother,self.z/fother)
    def __truediv__(self,other):
        try:
            fother = float(other)
        except:
            raise TypeError("Quaternions can only be divided by a real number.")
        return quaternion(self.w/fother,self.x/fother,self.y/fother,self.z/fother)
    def __pow__(self,other):
        """__pow__ will be prone to huge amounts of floating point issues. Be careful!"""
        try:
            fother = float(other)
        except:
            raise TypeError("Quaternions can only be raised to a real exponent.")
        #Next two lines are likely main sources of floating point errors
        theta = acos(self.w/abs(self))
        n_hat = quaternion(0,self.x,self.y,self.z)/abs(self)/sin(theta)
        return self.__norm__()**fother * (quaternion(cos(fother*theta),0,0,0) + n_hat*sin(fother*theta))
    def exp(self):
        """
        Returns the exponential of the quaternion.
        """
        vec = quaternion(0,self.x,self.y,self.z)
        norm_vec = abs(vec)
        return exp(self.w) * (quaternion(cos(norm_vec),0,0,0)+vec/norm_vec*sin(norm_vec))
    def log(self):
        """
        Returns the natural logarithm of the quaternion.
        """
        vec = quaternion(0,self.x,self.y,self.z)
        norm_vec = abs(vec)
        return quaternion(log(abs(self)),0,0,0) + vec/norm_vec *acos(self.w/abs(self))
