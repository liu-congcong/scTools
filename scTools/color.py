import ctypes

from .c import findLibrary


def generateColors(n):
    library = ctypes.cdll.LoadLibrary(findLibrary('color'))
    library.generateColors.argtypes = [ctypes.c_int]
    library.generateColors.restype = ctypes.POINTER(ctypes.c_double * 3)
    library.freeColors.argtypes = [ctypes.POINTER(ctypes.c_double * 3)]
    library.freeColors.restype = ctypes.c_int
    pointer = library.generateColors(n)
    assert pointer, 'Failed to generate colors.'
    colors = [(pointer[i][0], pointer[i][1], pointer[i][2]) for i in range(n)]
    library.freeColors(pointer)
    return colors
