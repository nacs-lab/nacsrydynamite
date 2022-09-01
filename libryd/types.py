from enum import Enum

class ParamType(Enum):
    Omega = 0
    Delta = 1
    C6 = 2
    Unknown = 3
    @staticmethod
    def fromStr(this_str):
        if "Omega".casefold() == this_str.casefold():
            return ParamType.Omega
        if "Delta".casefold() == this_str.casefold():
            return ParamType.Delta
        if "C6".casefold() == this_str.casefold():
            return ParamType.C6
        return ParamType.Unknown
    def __str__(self):
        return self.name
    def __repr__(self):
        return self.name
