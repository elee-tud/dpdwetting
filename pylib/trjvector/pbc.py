#!/home/elee/Programs/Python-3.6.10/bin/python


def shortest_vector_across_box(x, y, box):
    z=x-y
    for i in range(3):
        if z[i]>box[i]/2:
            z[i] -= box[i]
        elif z[i] < -box[i]/2.:
            z[i] += box[i]
    return z




