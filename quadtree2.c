
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <crtdbg.h>



const double epsilon_0 = 1.0e-3;

const double theta = 0.5;


/*a struct for storing particle information*/
typedef struct particle
{

    double x,y;
    
    double m;

    double v_x,v_y;

    double brightness;

    double f_x,f_y;

}particle;


/*struct to represent a quadtree node containing a particle as well as 4 subnodes*/
typedef struct node
{
    double x,y; //location of this node, as in where the four quadrants are split
    
    double half_dim; //half_dimension of this node

    double m; //mass of this node, for leaf node it is the mass of its particle; for stem node it is the sum of its particles

    double center_x, center_y; //location of the center of mass of this node

    double particle_x, particle_y; //location of the particle within the node (if it is leaf node)
    
    struct node *subnode_NW, *subnode_NE, *subnode_SW, *subnode_SE; //pointers to 4 subtree nodes

    bool stem; // wheather it is a leaf node or a stem node

}node;
/*note: the 4 subnodes are marked clock-wise, i.e. NW is 0, NE is 1, SE = 2 and SW = 3*/


/*function to construct a new node and initialize attributes */
void node_init(node** n, double x, double y, double half_dim) 
{

    
    (*n)->x = x;
    (*n)->y = y; //set location of quadrant meeting point

    (*n)->half_dim = half_dim; //quadrant "radius"


    (*n)->center_x = x;
    (*n)->center_y = y; //default center of mass set as the meeting point of 4 quadrants

    (*n)->stem = true; //stem node, not a leaf node
    (*n)->m = 0; // default no particles included

    (*n)->subnode_NE = NULL;
    (*n)->subnode_NW = NULL;
    (*n)->subnode_SE = NULL;
    (*n)->subnode_SW = NULL; 

    

    return;
}



double distance(node* stem, particle *p)
{


    double dist = sqrt(pow(stem->center_x - p->x,2) + pow(stem->center_y - p->y,2));
    return dist;

}


/*function calculating the force between a particle and the rest of the tree*/
void update_force (node** stem, particle* p, double G)
{

    if ((*stem)->stem == true) //if stem node, consider whether to use center of mass
    {
        double d = distance(*stem,p);
        if (2*(*stem)->half_dim/d > theta) //too close to be treated as single body
        {
            
                if ((*stem)->subnode_NE != NULL)
                {
                    update_force(&((*stem)->subnode_NE),p,G);
                }
                if ((*stem)->subnode_NW != NULL)
                {
                    update_force(&((*stem)->subnode_NW),p,G);
                }
                if ((*stem)->subnode_SE != NULL)
                {
                    update_force(&((*stem)->subnode_SE),p,G);
                }
                if ((*stem)->subnode_SW != NULL)
                {
                    update_force(&((*stem)->subnode_SW),p,G);
                }
            return;
        }

        else //can be considered as one body 
        {
            p->f_x = p->f_x - G*p->m*(*stem)->m/(pow(distance(*stem,p) + epsilon_0,3))*(p->x - (*stem)->center_x);
            p->f_y = p->f_y - G*p->m*(*stem)->m/(pow(distance(*stem,p) + epsilon_0,3))*(p->y - (*stem)->center_y);

            return;

        }
    }

    else //if it is a leaf node, simply calculate force between them
    {

        double r= sqrt(pow(p->x - (*stem)->particle_x,2) + pow(p->y - (*stem)->particle_y,2));

        p->f_x = p->f_x - G*p->m*(*stem)->m/(pow(r + epsilon_0,3))*(p->x - (*stem)->particle_x);
        p->f_y = p->f_y - G*p->m*(*stem)->m/(pow(r + epsilon_0,3))*(p->y - (*stem)->particle_y);
        
        return;
    }

    return;
}


/*function to update location of a particle after using its force*/
void update_location(particle *p, double timestep)
{

    double a_x = p->f_x/p->m;
    double a_y = p->f_y/p->m;

    p->v_x = p->v_x + a_x*timestep;
    p->v_y = p->v_y + a_y*timestep; //update velocity

    p->x = p->x + timestep*p->v_x;
    p->y = p->y + timestep*p->v_y; //update location
    
    return;

}


/*function to update mass and center of mass of a stem node when a new particle is added into a leaf node*/
void update_mass(node** stem, node* leaf)
{
    
    double m_sum = (*stem)->m + leaf->m;

    double x = ((*stem)->center_x*(*stem)->m + leaf->x*leaf->m)/m_sum;
    double y = ((*stem)->center_y*(*stem)->m + leaf->y*leaf->m)/m_sum;

    (*stem)->center_x = x;
    (*stem)->center_y = y; //update center of mass coordinates

    (*stem)->m = m_sum; //update stem node mass

 

    return;

}



void assign_particle(node** stem, particle *p)
{

    /*NW quadrant*/
    if (p->x < (*stem)->x && p->y > (*stem)->y)
    {

        /*NW leaf node is empty, insert particle*/
        if ((*stem)->subnode_NW == NULL)
        {
            node *NW = (node*)malloc(sizeof(node));
            NW->m = p->m;
            NW->particle_x = p->x;
            NW->particle_y = p->y;
            NW->stem = false;
            (*stem)->subnode_NW = NW;

            update_mass(&(*stem), ((*stem)->subnode_NW)); //update stem node center of mass
            
            

            return;
        }

        /*leaf node is not empty, must further subdivide*/
        else
        {
            particle *temp = (particle*)malloc(sizeof(particle));
            temp->m = (*stem)->subnode_NW->m;
            temp->x = (*stem)->subnode_NW->particle_x;
            temp->y = (*stem)->subnode_NW->particle_y; //store particle information in temporary particle
            
            node_init(&((*stem)->subnode_NW),(*stem)->x - 0.5*(*stem)->half_dim, (*stem)->y + 0.5*(*stem)->half_dim, 0.5*(*stem)->half_dim);
            //turn leaf node into new stem node

            assign_particle(&((*stem)->subnode_NW), temp);
            assign_particle(&((*stem)->subnode_NW), p); //recursively add new nodes into the new stem node

            //free(temp);

            return;


        }


    }


    /*NE quadrant*/
    if (p->x > (*stem)->x && p->y > (*stem)->y)
    {

        /*NE leaf node is empty, insert particle*/
        if ((*stem)->subnode_NE == NULL)
        {
            node *NE = (node*)malloc(sizeof(node));

            NE->m = p->m;
            NE->particle_x = p->x;
            NE->particle_y = p->y;
            NE->stem = false;

            (*stem)->subnode_NE = NE;
            

            update_mass(&(*stem), ((*stem)->subnode_NE)); //update stem node center of mass

           // free(NE);

            return;
        }

        /*leaf node is not empty, must further subdivide*/
        else
        {
            particle *temp = (particle*)malloc(sizeof(particle));

            temp->m = (*stem)->subnode_NE->m;
            temp->x = (*stem)->subnode_NE->particle_x;
            temp->y = (*stem)->subnode_NE->particle_y;
            
            node_init(&((*stem)->subnode_NE),(*stem)->x + 0.5*(*stem)->half_dim, (*stem)->y + 0.5*(*stem)->half_dim, 0.5*(*stem)->half_dim);
            //turn leaf node into new stem node

            assign_particle(&((*stem)->subnode_NE), temp);
            assign_particle(&((*stem)->subnode_NE), p); //recursively add new nodes into the new stem node

            //free(temp);
           
            return;


        }

    }

    /*SE quadrant*/
    if (p->x > (*stem)->x && p->y < (*stem)->y)
    {

        /*SE leaf node is empty, insert particle*/
        if ((*stem)->subnode_SE == NULL)
        {
            node* SE = (node*)malloc(sizeof(node));

            SE->m = p->m;
            SE->particle_x = p->x;
            SE->particle_y = p->y;
            SE->stem = false;

            (*stem)->subnode_SE = SE; 

            update_mass(&(*stem), ((*stem)->subnode_SE)); //update stem node center of mass
 
           // free(SE);

            return;
        }

        /*leaf node is not empty, must further subdivide*/
        else
        {
            particle *temp = (particle*)malloc(sizeof(particle));

            temp->m = (*stem)->subnode_SE->m;
            temp->x = (*stem)->subnode_SE->particle_x;
            temp->y = (*stem)->subnode_SE->particle_y;
            
            node_init(&((*stem)->subnode_SE),(*stem)->x + 0.5*(*stem)->half_dim, (*stem)->y - 0.5*(*stem)->half_dim, 0.5*(*stem)->half_dim);
            //turn leaf node into new stem node

            assign_particle(&((*stem)->subnode_SE), temp);
            assign_particle(&((*stem)->subnode_SE), p); //recursively add new nodes into the new stem node

            //free(temp);

            return;
        }

    }

    /*SW quadrant*/
    if (p->x < (*stem)->x && p->y < (*stem)->y)
    {

        /*SW leaf node is empty, insert particle*/
        if ((*stem)->subnode_SW == NULL)
        {
            node* SW = (node*)malloc(sizeof(node));
            SW->m = p->m;
            SW->particle_x = p->x;
            SW->particle_y = p->y;
            SW->stem = false;

            (*stem)->subnode_SW = SW;

            update_mass(&(*stem), ((*stem)->subnode_SW)); //update stem node center of mass

            
           // free(SW);

            return;
        }

        /*leaf node is not empty, must further subdivide*/
        else
        {
            particle *temp = (particle*)malloc(sizeof(particle));

            temp->m = (*stem)->subnode_SW->m;
            temp->x = (*stem)->subnode_SW->particle_x;
            temp->y = (*stem)->subnode_SW->particle_y; //store particle information in temporary particle
            
            node_init(&((*stem)->subnode_SW),(*stem)->x - 0.5*(*stem)->half_dim, (*stem)->y + 0.5*(*stem)->half_dim, 0.5*(*stem)->half_dim);
            //turn leaf node into new stem node

            assign_particle(&((*stem)->subnode_SW), temp);
            assign_particle(&((*stem)->subnode_SW), p); //recursively add new nodes into the new stem node
            
            //free(temp);

            return;


        }


    }
 
}


/*a function to free all the dynamically allocated nodes*/
void free_node(node** root)
{

    if((*root)->subnode_NE != NULL)
    {
        if ((*root)->subnode_NE->stem == true)
        {
            free_node(&(*root)->subnode_NE);
        }

        else free((*root)->subnode_NE);
    }

    if((*root)->subnode_NW != NULL)
    {
        if ((*root)->subnode_NW->stem == true)
        {
            free_node(&(*root)->subnode_NW);
        }

        else free((*root)->subnode_NW);
    }

    if((*root)->subnode_SE != NULL)
    {
        if ((*root)->subnode_SE->stem == true)
        {
            free_node(&(*root)->subnode_SE);
        }

        else free((*root)->subnode_SE);
    }

    if((*root)->subnode_SW != NULL)
    {
        if ((*root)->subnode_SW->stem == true)
        {
            free_node(&(*root)->subnode_SW);
        }

        else free((*root)->subnode_SW);
    }

    free(root);
    return;
}


int main(int argc, char* argv[])
{

     int N;
     int nstep;
     int graphics;
     double timestep;
    
     double G;

     char px[10], py[10], m[10],vx[10],vy[10],br[10];

     FILE* file;



    if(argc > 1)
    {

         N = (int)strtol(argv[1], NULL, 10);
         nstep =  (int)strtol(argv[3], NULL, 10);
         file = fopen(argv[2], "rb");
         timestep = strtod(argv[4], NULL);
         graphics = (int)strtol(argv[5], NULL, 10);       
         G = 100.0/N;

     }

      particle *p = (particle*)malloc(N*sizeof(particle));

    
    
    for (int i=0; i<N; i++)
    {
        fread(px,sizeof(px), 1, file);
        p[i].x = strtod(px, NULL);
        fread(py, sizeof(py), 1, file);
        p[i].y = strtod(py, NULL);
        fread(m, sizeof(m), 1, file);
        p[i].m = strtod(m, NULL);
        fread(vx, sizeof(vx), 1, file);
        p[i].v_x = strtod(vx, NULL);
        fread(vy, sizeof(vy), 1, file);
        p[i].v_y = strtod(vy, NULL);
        fread(br, sizeof(br), 1, file);
        p[i].brightness = strtod(br, NULL);

    }

    node* root;

    for (int t=0;t<nstep;t++)
    {
        
        root = (node*)malloc(sizeof(node));
        node_init(&root,0.5,0.5,0.5); //initiate root node

        for (int i=0;i<N;i++)
        {
            assign_particle(&root, &p[i]);

        } //add particles into root node


        for (int i=0;i<N;i++)
        {
            p[i].f_x = 0;
            p[i].f_y = 0; //initialize force to be 0
            
            update_force(&root, &(p[i]),G);
            

        }

        

        for(int i =0;i<N;i++)
        {

            update_location(&p[i], timestep);
            
        }

        
    }

    FILE *file2;
    file2 = fopen("result.gal", "wb");

    for (int i=0; i<N;i++)
    {
        fwrite(&(p[i].x), sizeof(p[i].x), 1, file2 );
        fwrite(&(p[i].y), sizeof(p[i].y), 1, file2 );
        fwrite(&(p[i].m), sizeof(p[i].m), 1, file2 );
        fwrite(&(p[i].v_x), sizeof(p[i].v_x), 1, file2 );
        fwrite(&(p[i].v_y), sizeof(p[i].v_y), 1, file2 );
        fwrite(&(p[i].brightness), sizeof(p[i].brightness), 1, file2 );
    }

    free_node(&root);

    free(p);
    
    
    return 0;
}






