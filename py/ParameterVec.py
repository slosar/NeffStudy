##
## Simple parameter vector class
##

#
# Simple class for dealing with Parameter.
# Parameter has a name, a value, an error and some bounds
# Names are also latex names.

class Parameter:
    def __init__(self,name,value):
        self.name=name
        self.value=value

    def sameParam(self,param2):
        return self.name==param2.name

    def setValue(self,val):
        self.value=val


class ParamList(list):

    def __init__(self, *args):
        list.__init__(self, *args)
        
    def __getslice__(self,i,j):
        return ParamList(list.__getslice__(self, i, j))
    def __add__(self,other):
        return ParamList(list.__add__(self,other))
    def __mul__(self,other):
        return ParamList(list.__mul__(self,other))

    def nameList(self):
        return [p.name for p in self]

    def valueList(self):
        return [p.value for p in self]

    def setValue(self,name,val):
        self[self.nameList().index(name)].setValue(val)

    def value(self,name):
        return self[self.nameList().index(name)].value
        
        
def DefaultParamList():
    pl=[Parameter('tau',0.07),
           Parameter('omegac',0.11987),Parameter('As',2.204e-9),Parameter('theta',1.0407781e-2),
           Parameter('Neff',3.046),Parameter('mnu',0.06),
           Parameter('omegab',0.022252),Parameter('ns',0.96475)]
    return ParamList(pl)
