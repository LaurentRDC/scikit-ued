# Copyright 2010 Torbjorn Bjorkman
# This file is part of cif2cell
#
# cif2cell is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cif2cell is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with cif2cell.  If not, see <http://www.gnu.org/licenses/>.
#
#******************************************************************************************
#  Description: Definitions of basic classes  and functions for CIF2Cell.
#
#  Author:      Torbjorn Bjorkman, torbjorn.bjorkman(at)aalto.fi
#  Affiliation: COMP, Aaalto University School of Science,
#               Department of Applied Physics, Espoo, Finland
#******************************************************************************************
from __future__ import division
from math import sqrt,acos
from elementdata import *
################################################################################################
# Miscellaneous
zero = 0.0
one = 1.0
two = 2.0
three = 3.0
four = 4.0
six = 6.0
third = 1/3
half = one/two
fourth = one/four
sixth = one/six
occepsilon = 0.000001
floatlist = [third, 2*third, half, fourth, one, zero, sqrt(2.0),sixth,5*sixth,sqrt(3.0),sqrt(3.0)/2]
angtobohr = 1.8897261
uperatogpercm = 1.6605388
uperautogpercm = 11.205871
#
codename = { 'abinit' : 'ABINIT',
             'ase' : 'ASE',
             'castep' : 'CASTEP',
             'cpmd' : 'CPMD',
             'cp2k' : 'CP2K',
             'cellgen' : 'cellgen',
             'crystal09' : 'CRYSTAL09',
             'elk' : 'Elk',
             'emto' : 'EMTO',
             'exciting' : 'exciting',
             'fleur' : 'Fleur',
             'ncol' : 'NCOL',
             'rspt' : 'RSPt',
             'siesta' : 'Siesta',
             'spacegroup' : 'spacegroup',
             'sprkkr' : 'SPR-KKR',
             'vasp' : 'VASP',
             'xband' : 'XBAND',
             'xyz' : 'xyz',
             'pwscf' : 'PWscf',
             'quantum-espresso' : 'Quantum Espresso'}
             
# Make a list of safe functions for use in safe_matheval. Thanks Lybniz developers!
safe_list = ['math','acos', 'asin', 'atan', 'atan2', 'ceil', 'cos', 'cosh', 'e', 'exp', 'fabs', 'floor', 'fmod', 'frexp', 'hypot', 'ldexp', 'log', 'log10', 'modf', 'pi', 'pow', 'radians', 'sin', 'sinh', 'sqrt', 'tan', 'tanh']
safe_dict = dict([ (k, locals().get(k, None)) for k in safe_list ])
# Add any needed builtins.
safe_dict['abs'] = abs

################################################################################################
# Exception classes
class SymmetryError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class PositionError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
class CellError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
class GeometryObjectError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
class SetupError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
################################################################################################
class GeometryObject:
    """
    Parent class for anything geometrical contains:
    compeps    : epsilon for determining when two floats are equal
    invcompeps : 1/compeps
    """
    def __init__(self,compeps=0.0002):
        self.compeps = compeps
        self.invcompeps = 1./self.compeps

class CellFloat(float, GeometryObject):
    """
    Hashable float with comparison inherited from GeometryObject.
    """
    def __init__(self,t):
        float.__init__(self)
        GeometryObject.__init__(self)
    def __hash__(self):
        return int(self*self.invcompeps)
    def __eq__(self, other):
        if abs(self-other) > self.compeps:
            return False
        else:
            return True
        
class Charge(float):
    """
    Class for representing the charge state/oxidation number of an atom (ion).
    It is just an integer, except for a modified routine for turning into a string,
    viz. a two plus ion gets the string representation '2+'
    """
    def __init__(self,i):
        float.__init__(self)
    def __str__(self):
        if abs(self-int(self)) < 0.0001:
            if int(self) == 0:
                return '0'
            elif self > 0:
                return str(abs(int(self)))+'+'
            elif self < 0:
                return str(abs(int(self)))+'-'
        else:
            if self > 0:
                return str(abs(self))+"+"
            else:
                return str(abs(self))+"-"

class Vector(list,GeometryObject):
    """
    Floating point vector describing a lattice point. Assumed to have length three.
    * Supports testing for equality using the self.compeps parameter inherited
    from the parent class GeometryObject
    * Supports testing for < using the euclidean norm (length) of the vector
    * Supports addition with another Vector, in which case the
      vectors are added component-wise.
    More methods:
    length            : returns the euclidean norm of the vector
    transform(matrix) : returns matrix*vector
    improveprecision  : identify some conspicuous numbers and improve precision
    """
    def __init__(self, vec, *args, **kwargs):
        GeometryObject.__init__(self, *args, **kwargs)
        list.__init__(self, [float(v) for v in vec])
    def __hash__(self):
        t = self.invcompeps*self[0]+1.1*self.invcompeps*self[1]+1.2*self.invcompeps*self[2]
        return int(round(t))
    def __eq__(self,other):
        for i in range(3):
            if abs(self[i]-other[i]) > self.compeps:
                return False
        return True
    def __lt__(self, other):
        sl = self[0]**2+self[1]**2+self[2]**2
        ol = other[0]**2+other[1]**2+other[2]**2
        return sl < ol
    # Addition of two vectors
    def __add__(self, other):
        t = Vector([self[i]+other[i] for i in range(3)])
        return t
    # Subtraction of two vectors
    def __sub__(self, other):
        t = Vector([self[i]-other[i] for i in range(3)])
        return t
    def __neg__(self):
        return Vector([-t for t in self])
    def __str__(self):
        s = ""
        for e in self:
            if type(e) == type(1):
                s += "%2i "%e
            else:
                s+= "%19.15f "%e
        return s
    # Length of the vectors
    def length(self):
        return sqrt(self[0]**2+self[1]**2+self[2]**2)
    # Multiplication by scalar
    def scalmult(self, a):
        t = []
        for i in range(3):
            t.append(self[i]*a)
        return Vector(t)
    # dot product
    def dot(self,a):
        t = 0.0
        for i in range(3):
            t += self[i]*a[i]
        return t
    # triple product
    def triple(self,a,b):
        t = [a,b,c]
        return det3(t)
    # coordinate transformation
    def transform(self, matrix):
        t = Vector(mvmult3(matrix, self))
        return t
    def improveprecision(self):
        for i in range(3):
            self[i] = improveprecision(self[i],self.compeps)
        return self
    # Angle between this vector and another vector
    def angle(self, other):
        return acos(self.dot(other)/(self.length() * other.length()))

class LatticeVector(Vector):
    """
    Vector of length three that maps back things into the cell
    """
    def __init__(self, vec, interval=(0.0, 1.0), *args, **kwargs):
        Vector.__init__(self, vec, *args, **kwargs)
        # Interval we wish to use for the coordinates.
        # In practice either [0,1] or [-.5, 0.5]
        self.interval = interval
        self.improveprecision()
        self.intocell()
    # Addition of two vectors, putting the result back
    # into the cell if necessary
    def __add__(self, other):
        if self.interval[0] != other.interval[0] or self.interval[1] != other.interval[1]:
            raise GeometryObjectError("LatticeVectors must have the same definition intervals to be added.")
        t = LatticeVector([self[i]+other[i] for i in range(3)])
        t.intocell()
        return t
    # Change interval
    def change_interval(self, interval):
        self.interval = interval
        t = LatticeVector([0,0,0],interval=interval)
        self = self + t
    # coordinate transformation
    def transform(self, matrix):
        t = LatticeVector(mvmult3(matrix, self))
        t.intocell()
        return t
    # Put the vector components into the cell interval defined by self.interval
    def intocell(self):
        for i in range(3):
            while self[i] < self.interval[0]:
                self[i] += 1.0
            while self[i] >= self.interval[1]-self.compeps:
                self[i] -= 1.0

class LatticeMatrix(GeometryObject, list):
    """
    Three by three matrix
    """
    def __init__(self,mat):
        GeometryObject.__init__(self)
        t = []
        for vec in mat:
            t.append(Vector(vec))
        list.__init__(self, t)
    # no idea whether this is a clever choice of hash function...
    def __hash__(self):
        t = self.invcompeps*(self[0][0] + self[0][1] + self[0][2] +\
                             self[0][0] + self[0][1] + self[0][2] +\
                             self[0][0] + self[0][1] + self[0][2])
        return int(t)
    def __str__(self):
        matstr = ""
        for l in self:
            matstr += str(l)+"\n"
        return matstr
    def __eq__(self,other):
        for i in range(3):
            for j in range(3):
                if abs(self[i][j]-other[i][j]) > self.compeps:
                    return False
        return True
    # coordinate transformation
    def transform(self, matrix):
        return LatticeMatrix(mmmult3(matrix, self))
    # transpose
    def transpose(self):
        t = [[self[0][0], self[1][0], self[2][0]],
             [self[0][1], self[1][1], self[2][1]],
             [self[0][2], self[1][2], self[2][2]]]
        return LatticeMatrix(t)
    # Round floats to recognized numbers
    def improveprecision(self):
        t = []
        for i in range(3):
            t.append([])
            for j in range(3):
                t[i].append(improveprecision(self[i][j],self.compeps))
        return LatticeMatrix(t)
                            
class AtomSite(GeometryObject):
    """
    Class for describing an atomic site.

    Contains data:
    position  : a vector that gives the position
    species   : a dictionary with element-occupancy pairs (e.g. {Fe : 0.2, Co : 0.8})
    label     : any label
    charges   : a dictionary with the charge states (oxidation numbers) of the different species
    index     : any integer

    Functions:
        __eq__    : compare equality
        __str__   : one line with species and position info
        copy      : return a new instance with the same data
        distance  : distance to some other atom
        spcstring : species string ('Mn', 'La/Sr' ...)
        alloy     : true if there are more than one species occupying the site
        
    """
    def __init__(self,position=None,species=None,label="",charges=None,index=None):
        GeometryObject.__init__(self)
        if position != None:
            self.position = LatticeVector(position)
        else:
            self.position = None
        if species != None:
            self.species = species
        else:
            self.species = {}
        if charges != None:
            self.charges = charges
        else:
            if self.species != None:
                self.charges = {}
                for k in self.species.keys():
                    self.charges[k] = Charge(0)
            else:
                self.charges = {}
        self.label = label
        ## self.charge = Charge(charge)
        self.index = index
    def __hash__(self):
        return hash(hash(self.position)+hash(''.join(sorted(self.species.keys())))+hash(sum(self.species.values())))
    def __eq__(self,other):
        return self.position == other.position and self.species == other.species
    # Species string, sorted by atomic weight.
    def spcstring(self):
        t = []
        for k in self.species:
            t.append(k)
        t.sort(key = lambda x: ElementData().elementweight[x], reverse=True)
        tmp = ""
        for k in t:
            tmp += k+"/"
        tmp = tmp.rstrip("/")
        return tmp
    # Is there more than one species on this site?
    def alloy(self):
        occsum = sum([v for k,v in self.species.iteritems()])
        return len(self.species) > 1 or abs(occsum-1) > self.compeps
    # print site data in some informative way
    def __str__(self):
        # Element symbol
        tmp = self.spcstring().ljust(8)
        # Position
        tmp += " %19.15f %19.15f %19.15f   "%(self.position[0],self.position[1],self.position[2])
        # occupancy
        for k,v in self.species.iteritems():
            tmp += str(v)+"/"
        tmp = tmp.rstrip("/")
        return tmp
    def CIradius(self,size="max",covalent=False):
        """
        Return maximal/minimal Covalent/Ionic radius of the site.
        'size' controls whether the maximal or minimal radius is returned
        'covalent' will enforce the covalent radius.
        """
        t = []
        if covalent:
            for sp in self.species.keys():
                try:
                    t.append(ElementData().CovalentRadius2[sp])
                except:
                    pass
        else:
            for sp,ch in self.charges.iteritems():
                try:
                    t.append(ElementData().IonicRadius[sp+str(ch)])
                except:
                    try:
                        t.append(ElementData().CovalentRadius2[sp])
                    except:
                        pass
        try:
            if size == "min":
                return min(t)
            else:
                return max(t)
        except:
            return None
    #
    # The distance to another atom site
    def distance(self,other):
        v = self.position-other.position
        return v.length()

class SymmetryOperation(GeometryObject):
    """
    Class describing a symmetry operation, with a rotation matrix and a translation.
    """
    def __init__(self, eqsite=None):
        GeometryObject.__init__(self)
        self.eqsite = eqsite
        if self.eqsite != None:
            self.rotation = self.rotmat()
            self.translation = LatticeVector(self.transvec())
        else:
            self.rotation = None
            self.translation = None
    def __hash__(self):
        return hash(self.rotation)+hash(self.translation)
    # This way of printing was useful for outputting to CASTEP.
    def __str__(self):
        return str(self.rotation)+str(self.translation)+"\n"
    # Two symmetry operations are equal if rotation matrices and translation vector
    # differ by at most compeps
    def __eq__(self, other):
        eq = True
        for i in range(3):
            for j in range(3):
                eq = eq and abs(self.rotation[i][j] - other.rotation[i][j]) < self.compeps
            eq = eq and self.translation == other.translation
        return eq
    # Comparison between operations made by comparing lengths of translation vectors,
    # whether the rotation is diagonal and the identity is always less than anything else.
    # That way we only need to sort a list of operations to get identity first (and a reasonably
    # intuitive list order).
    def __lt__(self, other):
        if self.translation < other.translation:
            return True
        if other.translation < self.translation:
            return False
        if self.diagonal():
            # diagonal matrices "smaller"
            if not other.diagonal():
                return True
            # identity is "smallest"
            if self.rotation[0][0] == self.rotation[1][1] == self.rotation[2][2] == 1:
                return True
            return False
        else:
            return False
        return self.translation < other.translation
    # Return a rotation matrix from "x,y,z" representation of a symmetry operation
    # !!!With respect to cartesian axes!!!
    def rotmat(self):
        mat = [[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]
        for j in range(len(self.eqsite)):
            xyz = self.eqsite[j].replace('+',' +').replace('-',' -').split()
            for i in xyz:
                if i.strip("+-") == 'x':
                    mat[0][j] = float(i.strip('x')+"1")
                elif i.strip("+-") == 'y':
                    mat[1][j] = float(i.strip('y')+"1")
                elif i.strip("+-") == 'z':
                    mat[2][j] = float(i.strip('z')+"1")            
        return LatticeMatrix(mat)
    # Return a translation vector from "x,y,z" representation of a symmetry operation
    def transvec(self):
        vec = []
        for i in range(3):
            vec.append(0.0)
        for j in range(len(self.eqsite)):
            xyz = self.eqsite[j].replace('+',' +').replace('-',' -').split()
            for i in xyz:
                if i.strip("+-xyz") != "":
                    vec[j] = safe_matheval(i)
        return LatticeVector(vec)
    # True if the operation is diagonal
    def diagonal(self):
        if abs(self.rotation[0][1]) < self.compeps and \
           abs(self.rotation[0][2]) < self.compeps and \
           abs(self.rotation[1][0]) < self.compeps and \
           abs(self.rotation[1][2]) < self.compeps and \
           abs(self.rotation[2][0]) < self.compeps and \
           abs(self.rotation[2][1]) < self.compeps:
            return True
        else:
            return False
    # Operate on a vector and return the result
    def operate(self,vector):
        t = Vector(mvmult3(self.rotation, vector)) + self.translation
        return t
    def improveprecision(self):
        self.rotation = self.rotation.improveprecision()
        self.translation = self.translation.improveprecision()

################################################################################################
# Dictionaries
settingname = { 'P' : 'primitive',
                'A' : 'base-centered',
                'B' : 'base-centered',
                'C' : 'base-centered',
                'I' : 'body-centered',
                'F' : 'face-centered',
                'R' : 'rhombohedral' }
        
################################################################################################
# Functions
# Evaluate expr safely, i.e. only allow execution of mathematical functions
def safe_matheval(expr):
    sexpr = expr.replace(" ","")
    return eval(sexpr,{"__builtins__":None},safe_dict)

def removeerror(string):
    # Remove error estimates at the end of a number (as in 3.28(5))
    splitstr=string.split('(')
    return splitstr[0]

# Guess the "true" values of some conspicuous numbers
def improveprecision(x,eps):
    for f in floatlist:
        if abs(abs(x)-f) <= eps:                
            # 0
            return copysign(f,x)
    # if no match found, return x
    return x

def latvectadd(a,b):
    t = []
    for i in range(3):
        t.append(a[i]+b[i])
        if abs(t[i]) >= 1-occepsilon:
            t[i] = t[i] - copysign(1,t[i])
        t[i] = improveprecision(t[i],occepsilon)
    return t

def putincell(coords,coordepsilon):
    # Put coordinates in the interval 0 <= x < 1
    for i in range(3):
        # first make the coordinate positive
        while coords[i] < 0:
            coords[i] = coords[i] + 1
        # then put it in the primitive cell
        while coords[i] > 1-coordepsilon:
            coords[i] = coords[i] - 1
    return coords

# Determinant of 3x3 dimensional matrix
def det3(m):
    a = m[1][1]*m[2][2]-m[1][2]*m[2][1]
    b = m[1][2]*m[2][0]-m[1][0]*m[2][2]
    c = m[1][0]*m[2][1]-m[1][1]*m[2][0]
    return m[0][0]*a + m[0][1]*b + m[0][2]*c

# Inverse of 3x3 dimensional matrix
def minv3(m):
    di = 1/det3(m)
    w = [[(m[1][1]*m[2][2]-m[1][2]*m[2][1])*di, (m[0][2]*m[2][1]-m[0][1]*m[2][2])*di, (m[0][1]*m[1][2]-m[0][2]*m[1][1])*di],
         [(m[1][2]*m[2][0]-m[1][0]*m[2][2])*di, (m[0][0]*m[2][2]-m[0][2]*m[2][0])*di, (m[0][2]*m[1][0]-m[0][0]*m[1][2])*di],
         [(m[1][0]*m[2][1]-m[1][1]*m[2][0])*di, (m[0][1]*m[2][0]-m[0][0]*m[2][1])*di, (m[0][0]*m[1][1]-m[0][1]*m[1][0])*di]]
    return w

# matrix-vector multiplication
def mvmult3(mat,vec):
    w = [0.,0.,0.]
    for i in range(3):
        t = 0
        for j in range(3):
            t = t + mat[j][i]*vec[j]
        w[i] = t
    return w

# more efficient, but goes the other way...
## def mvmult3(mat,vec):
##     w = [0.,0.,0.]
##     t = 0
##     for i in range(3):
##         r = mat[i]
##         for j in range(3):
##             t += r[j]*vec[j]
##         w[i],t = t,0
##     return w

# matrix-matrix multiplication
def mmmult3(m1,m2):
    w = []
    for i in range(3):
        w.append([])
        for j in range(3):
            t = 0
            for k in range(3):
                t += m1[i][k]*m2[k][j]
            w[i].append(t)
    return w

def crystal_system(spacegroupnr):
    # Determine crystal system
    if 0 < spacegroupnr <= 2:
        return "triclinic"
    elif 2 < spacegroupnr <=15:
        return "monoclinic"
    elif 15 < spacegroupnr <= 74:
        return "orthorhombic"
    elif 74 < spacegroupnr <= 142:
        return "tetragonal"
    elif 142 < spacegroupnr <= 167:
        return "trigonal"
    elif 167 < spacegroupnr <= 194:
        return "hexagonal"
    elif 194 < spacegroupnr <= 230:
        return "cubic"
    else:
        return "unknown"

# Return x with the sign of y
def copysign(x, y):
    if y >= 0:
        return x
    else:
        return -x

# Get rid of newline characters from a string, portably.
# replace allows to specify a string to replace the newline
# character with.
def deletenewline(string, replace=""):
    tmp = string
    # Windows...
    if "\r\n" in string:
        tmp = tmp.replace("\r\n",replace)
    # ...Mac...
    if "\r" in string:
        tmp = tmp.replace("\r",replace)
    # ...and Unix.
    if "\n" in string:
        tmp = tmp.replace("\n",replace)
    return tmp
