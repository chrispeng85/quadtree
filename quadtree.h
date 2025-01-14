#ifndef QUADTREE_H
#define QUADTREE_H

#include <stdbool.h>

extern const double epsilon_0;
extern const double theta;

typedef struct particle {

    double x,y;
    double m;
    double v_x, v_y;
    double brighteness;
    double f_x, f_y;


} particle;


typedef struct node {

    double x, y;
    double half_dim;
    double m;
    double center_x, center_y;
    double particle_x, particle_y;
    struct node *subnode_NW, *subnode_NE, *subnode_SW, *subnode_SE;
    bool stem;

} node;


void node_init(node** n, double x, double y, double haflf_dim);
double distance(node* stem, particle *p);
void update_force(node** stem, particle* p, double G);
void update_location(particle *p, double *timestep);
void update_mass(node** stem, node* leaf);
void assign_particle(node** stem, particle *p);
void free_node(node* root);


#endif