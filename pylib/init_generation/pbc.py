import sys
import numpy as np

def get_particle_in_box(r, box):
    return np.array([get_particle_vector_in_box(r[0], box[0]), get_particle_vector_in_box(r[1], box[1]), get_particle_vector_in_box(r[2], box[2])])


def get_particle_vector_in_box(x, boxx):
    newx=x
    if newx<0:
        while newx<0:
            newx += boxx
    elif newx>=boxx:
        while newx>=boxx:
            newx -= boxx
    if newx<0:
        print(newx)
    return newx


def get_shortest_image_vector(r, box):
    return np.array([get_shortest_image_dist(r[0], box[0]), get_shortest_image_dist(r[1], box[1]),get_shortest_image_dist(r[2], box[2])])


def get_shortest_image_dist(x, boxx):
    newx=x
    halfx=boxx/2
    if newx>=halfx:
        newx-=boxx
    elif newx<-halfx:
        newx+=boxx
    return newx


