import ctypes
from datetime import datetime

from .c import findLibrary


def generateColors(n):
    library = ctypes.cdll.LoadLibrary(findLibrary('color'))
    library.generateColors.argtypes = [ctypes.c_int]
    library.generateColors.restype = ctypes.POINTER(ctypes.c_double * 3)
    library.freeColors.argtypes = [ctypes.POINTER(ctypes.c_double * 3)]
    library.freeColors.restype = ctypes.c_int
    pointer = library.generateColors(n)
    if pointer is None:
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Error: Failed to generate colors.', file = stderr, flush = True)
        exit(1)
    colors = [(pointer[i][0], pointer[i][1], pointer[i][2]) for i in range(n)]
    library.freeColors(pointer)
    return colors