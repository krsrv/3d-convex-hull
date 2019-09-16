/*
References :
1. Computational Geometry - Algorithms and Applications : Berg
2. http://www.cs.sfu.ca/~binay/813.2011/DCEL.pdf
*/

#include <list>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <set>
#include <unordered_set>
#include <iostream>

enum shape { NONE, LINEAR, PLANAR, TETRA };
struct vertex; struct face; struct halfEdge;

struct vertex{
  float x, y, z;
  halfEdge* inc;
};

struct face{
  halfEdge* boundary;
};

struct halfEdge{
  vertex* origin;
  halfEdge *twin, *next, *prev;
  face* incident;
};

std::map<vertex*, std::unordered_set<face*> > vertexConflict;
std::map<face*, std::unordered_set<vertex*> > faceConflict;

std::unordered_set<face*> CH;
vertex* interior;

bool verbose = true;

std::string coordinate(vertex *v){
  return "V(" + std::to_string((int)v->x) + ", " + std::to_string((int)v->y) + ", " + std::to_string((int)v->z) + ")";
}

std::string faces(face *f){
  vertex* fv[3] = {f->boundary->origin, f->boundary->next->origin, f->boundary->next->next->origin};
  return "F[" + coordinate(fv[0]) + ", " + coordinate(fv[1]) + ", " + coordinate(fv[2]) + "]";
}

std::string edges(halfEdge* e){
  vertex* fv[2] = {e->origin, e->twin->origin};
  return "E[" + coordinate(fv[0]) + ", " + coordinate(fv[1]) + "]";
}

void detail(){
  std::set<vertex*> set;
  std::cout << CH.size() << " Faces:\n";
  for(std::unordered_set<face*>::iterator it = CH.begin(); it != CH.end(); ++it){
    std::cout << faces(*it) << std::endl;
    std::cout << "Edges:\n\t" << edges((*it)->boundary) << "\n";
    set.insert((*it)->boundary->origin);
    for(halfEdge* i = (*it)->boundary->next; i != (*it)->boundary; i = i->next){
      std::cout << "\t" << edges(i) << "\n";
      set.insert(i->origin);
    }
  }
  std::cout << set.size() << std::endl;
  for(std::set<vertex*>::iterator it = set.begin(); it != set.end(); ++it)
    std::cout << (*it)->x << " " << (*it)->y << " " << (*it)->z << std::endl;
}

int makeTetrahedron(std::vector<vertex*> &input, int length){
  if(length == 1) return NONE;
  else if(length == 2) return LINEAR;
  vertex *v[4];
  v[0] = input[0];
  v[1] = input[1];
  v[2] = v[3] = input[0] = input[1] = nullptr;
  
  int i;
  vertex b1 = {
    .x = v[1]->x - v[0]->x,
    .y = v[1]->y - v[0]->y,
    .z = v[1]->z - v[0]->z,
  };
  
  for(i = 2; i < length; i++){
    vertex b2 = {
      .x = input[i]->x - v[0]->x,
      .y = input[i]->y - v[0]->y,
      .z = input[i]->z - v[0]->z
    };
    if(b2.y*b1.z-b2.z*b1.y == 0 && b2.z*b1.x-b2.x*b1.z == 0 && b2.x*b1.y-b2.y*b1.x == 0)  continue;
    v[2] = input[i];
    input[i] = nullptr;
    break;
  }
  if(v[2] == nullptr) return LINEAR;
 
  float normal[3] = {
                      (v[2]->y - v[0]->y)*(v[1]->z - v[0]->z) - (v[2]->z - v[0]->z)*(v[1]->y - v[0]->y),
                      (v[2]->z - v[0]->z)*(v[1]->x - v[0]->x) - (v[2]->x - v[0]->x)*(v[1]->z - v[0]->z),
                      (v[2]->x - v[0]->x)*(v[1]->y - v[0]->y) - (v[2]->y - v[0]->y)*(v[1]->x - v[0]->x),
                    };
  for(i = i+1; i < length; i++){
    vertex b2 = {
      .x = input[i]->x - v[0]->x,
      .y = input[i]->y - v[0]->y,
      .z = input[i]->z - v[0]->z
    };
    if(b2.x * normal[0] + b2.y * normal[1] + b2.z * normal[2] == 0) continue;
    v[3] = input[i];
    input[i] = nullptr;
    break;
  }
  if(v[3] == nullptr) return PLANAR;

  // Initialise the tetrahedron
  float norm[3] = {
    (v[2]->y - v[0]->y)*(v[1]->z - v[0]->z) - (v[2]->z - v[0]->z)*(v[1]->y - v[0]->y),
    (v[2]->z - v[0]->z)*(v[1]->x - v[0]->x) - (v[2]->x - v[0]->x)*(v[1]->z - v[0]->z),
    (v[2]->x - v[0]->x)*(v[1]->y - v[0]->y) - (v[2]->y - v[0]->y)*(v[1]->x - v[0]->x),
  };
  vertex b = {
    .x = v[3]->x - v[0]->x,
    .y = v[3]->y - v[0]->y,
    .z = v[3]->z - v[0]->z
  };
  
  if(b.x * norm[0] + b.y * norm[1] + b.z * norm[2] < 0){
    vertex* &buff = v[0];
    v[0] = v[1];
    v[1] = buff;
  }

  if(verbose) std::cout << "Initialising tetrahedron\n";
  halfEdge *e01 = new halfEdge, *e02 = new halfEdge, *e03 = new halfEdge, *e12 = new halfEdge, *e13 = new halfEdge, *e23 = new halfEdge, *e10 = new halfEdge, *e20 = new halfEdge, *e30 = new halfEdge, *e21 = new halfEdge, *e31 = new halfEdge, *e32 = new halfEdge;
  e01->origin = e02->origin = e03->origin = v[0];
  e10->origin = e12->origin = e13->origin = v[1];
  e20->origin = e21->origin = e23->origin = v[2];
  e30->origin = e31->origin = e32->origin = v[3];
  e01->twin = e10;  e10->twin = e01;
  e02->twin = e20;  e20->twin = e02;
  e03->twin = e30;  e30->twin = e03;
  e12->twin = e21;  e21->twin = e12;
  e13->twin = e31;  e31->twin = e13;
  e23->twin = e32;  e32->twin = e23;
  e20->next = e01;  e01->next = e12;  e12->next = e20;
  e02->next = e23;  e23->next = e30;  e30->next = e02;
  e10->next = e03;  e03->next = e31;  e31->next = e10;
  e21->next = e13;  e13->next = e32;  e32->next = e21;
  e01->prev = e20;  e20->prev = e12;  e12->prev = e01;
  e23->prev = e02;  e02->prev = e30;  e30->prev = e23;
  e03->prev = e10;  e10->prev = e31;  e31->prev = e03;
  e13->prev = e21;  e21->prev = e32;  e32->prev = e13;
  face *f[4];
  f[0] = new face; f[1] = new face; f[2] = new face; f[3] = new face;
  f[0]->boundary = e01;
  f[1]->boundary = e02;
  f[2]->boundary = e03;
  f[3]->boundary = e13;
  e01->incident = e12->incident = e20->incident = f[0];
  e02->incident = e23->incident = e30->incident = f[1];
  e03->incident = e31->incident = e10->incident = f[2];
  e13->incident = e32->incident = e21->incident = f[3];
  v[0]->inc = e01;
  v[1]->inc = e12;
  v[2]->inc = e23;
  v[3]->inc = e31;
  CH.insert(f[0]);
  CH.insert(f[1]);
  CH.insert(f[2]);
  CH.insert(f[3]);

  interior = new vertex;
  interior->x = 0.25*(v[0]->x + v[1]->x + v[2]->x + v[3]->x);
  interior->y = 0.25*(v[0]->y + v[1]->y + v[2]->y + v[3]->y);
  interior->z = 0.25*(v[0]->z + v[1]->z + v[2]->z + v[3]->z);

  v[0] = v[1] = v[2] = v[3] = nullptr;

  return TETRA;
}

// Check whether two faces are coplanar using their normals
bool coplanar(face* f, face* g){
  vertex* fv[3] = {f->boundary->origin, f->boundary->next->origin, f->boundary->next->next->origin};
  float normal1[3] = {
                      (fv[2]->y - fv[0]->y)*(fv[1]->z - fv[0]->z) - (fv[2]->z - fv[0]->z)*(fv[1]->y - fv[0]->y),
                      (fv[2]->z - fv[0]->z)*(fv[1]->x - fv[0]->x) - (fv[2]->x - fv[0]->x)*(fv[1]->z - fv[0]->z),
                      (fv[2]->x - fv[0]->x)*(fv[1]->y - fv[0]->y) - (fv[2]->y - fv[0]->y)*(fv[1]->x - fv[0]->x),
                    };
  vertex* gv[3] = {g->boundary->origin, g->boundary->next->origin, g->boundary->next->next->origin};
  float normal2[3] = {
                      (gv[2]->y - gv[0]->y)*(gv[1]->z - gv[0]->z) - (gv[2]->z - gv[0]->z)*(gv[1]->y - gv[0]->y),
                      (gv[2]->z - gv[0]->z)*(gv[1]->x - gv[0]->x) - (gv[2]->x - gv[0]->x)*(gv[1]->z - gv[0]->z),
                      (gv[2]->x - gv[0]->x)*(gv[1]->y - gv[0]->y) - (gv[2]->y - gv[0]->y)*(gv[1]->x - gv[0]->x),
                    };
  
  fv[0] = fv[1] = fv[2] = gv[0] = gv[1] = gv[2] = nullptr;
//  std::cout << faces(f) << std::endl << faces(g) << std::endl;
//  std::cout << normal1[0] << ", " << normal1[1] << ", " << normal1[2] << "\n";
//  std::cout << normal2[0] << ", " << normal2[1] << ", " << normal2[2] << "\n";
  if(normal1[0] == normal2[0] && normal1[1] == normal2[1] && normal1[2] == normal2[2])  return true;
  return false;
}

// Check whether the face is visible from the given vertex
bool isVisible(face* f, vertex* v){
  vertex* fv[3] = {f->boundary->origin, f->boundary->next->origin, f->boundary->next->next->origin};
  float normal[3] = {
                      (fv[2]->y - fv[0]->y)*(fv[1]->z - fv[0]->z) - (fv[2]->z - fv[0]->z)*(fv[1]->y - fv[0]->y),
                      (fv[2]->z - fv[0]->z)*(fv[1]->x - fv[0]->x) - (fv[2]->x - fv[0]->x)*(fv[1]->z - fv[0]->z),
                      (fv[2]->x - fv[0]->x)*(fv[1]->y - fv[0]->y) - (fv[2]->y - fv[0]->y)*(fv[1]->x - fv[0]->x),
                    };
  vertex iV = {
    .x = interior->x - fv[0]->x,
    .y = interior->y - fv[0]->y,
    .z = interior->z - fv[0]->z
  };
  vertex eV = {
    .x = v->x - fv[0]->x,
    .y = v->y - fv[0]->y,
    .z = v->z - fv[0]->z
  };
 
  fv[0] = fv[1] = fv[2] = nullptr;

  float norm[2] = {eV.x*normal[0]+eV.y*normal[1]+eV.z*normal[2], iV.x*normal[0]+iV.y*normal[1]+iV.z*normal[2]};
  if(norm[0] * norm[1] >= 0) return false;
  else                      return true;
}

// Initialise the conflict graph
void makeConflict(std::vector<vertex*> &v, int length){
  for(int i = 2; i < length; i++){
    if(v[i] == nullptr) continue;
    for(std::unordered_set<face*>::iterator j = CH.begin(); j != CH.end(); ++j)
      if(isVisible(*j, v[i])){
        if(verbose){
          std::cout << "Conflict found : " << coordinate(v[i]) << " and " << faces(*j) << "\n";
        }
        faceConflict[*j].insert(v[i]);
        vertexConflict[v[i]].insert(*j);
      }
  }
}

void deleteAll(std::vector<face*> &fList, vertex* v){
  if(verbose) std::cout << "Deleting edges of " << fList.size() << "\n";
  for(std::vector<face*>::iterator i = fList.begin(); i != fList.end(); ++i){
    halfEdge* start = (*i)->boundary;
    for(halfEdge* it = start->next; it != start;){
      halfEdge *jt = it->next;
      std::cout << "Deleting " << std::flush << coordinate(it->origin) << std::flush << "\n";
      delete it;
      it = jt;
    }
    delete start;

    if(faceConflict.find(*i) != faceConflict.end()){
      for(std::unordered_set<vertex*>::iterator j = faceConflict[*i].begin(); j != faceConflict[*i].end(); ++j){
            std::unordered_set<face*>::iterator fa = vertexConflict[*j].find(*i);
            if(fa != vertexConflict[*j].end())  vertexConflict[*j].erase(fa);
      }
    }
    std::map<face*, std::unordered_set<vertex*> >::iterator fa = faceConflict.find(*i);
    if(fa != faceConflict.end())  faceConflict.erase(fa);
    //fList->second.end();
  }
  
  std::map<vertex*, std::unordered_set<face*> >::iterator it = vertexConflict.find(v);
  if(it != vertexConflict.end())  vertexConflict.erase(it);
  
}

// Obtain the horizon list for the given vertex and polyhedron using the list of conflicting faces. Simultaneously delete the conflicitng faces
void obtainHorizon(vertex* v, std::unordered_set<face*> &fList, std::vector<halfEdge*> &horizon){
  halfEdge *start;
  bool BREAK;
  for(std::unordered_set<face*>::iterator it = fList.begin(); it != fList.end(); ++it){
    start = (*it)->boundary;
    if(!isVisible(start->twin->incident, v)) break;
    for(start = start->next; start != (*it)->boundary ; start = start->next)
      if(!isVisible(start->twin->incident, v)){ 
        BREAK = true;
        break;
      }
    if(BREAK) break;
  }
  horizon.push_back(start->twin);
  if(verbose) std::cout << "Added " << edges(start->twin) << " to Horizon List\n";

  halfEdge* edge = start->twin->next;
  while(true){
    while(true){
      bool x = isVisible(edge->incident, v), y = isVisible(edge->twin->incident, v);
      if(x != y)  break;
      edge = edge->twin->next;
    }
    if(edge == start->twin) break;
    horizon.push_back(edge);
    if(verbose) std::cout << "Added " << edges(edge) << " to Horizon List\n";
    edge = edge->next;
  }
  start = edge = nullptr;
}

// Create the 3D Convex Hull
void TetraHull(std::vector<vertex*> &v, int length){
  detail();
  for(int i = 2; i < length; i++){
    if(v[i] == nullptr) continue;
    std::map<vertex*, std::unordered_set<face*> >::iterator vi = vertexConflict.find(v[i]);
    if(vi == vertexConflict.end()) continue;
    if(vi->second.size() == 0){
      vertexConflict.erase(vi);
      continue;
    }
    
    if(verbose) std::cout << "Conflicting facets exist for Vertex " << coordinate(v[i]) << "\n";
    // Conflicting facets exist. Obtain the list through the conflict graph
    std::unordered_set<face*> facetList = vi->second;
    for(std::unordered_set<face*>::iterator j = facetList.begin(); j != facetList.end(); ++j){
      // Remove the conflicting facet from the current Convex Hull
      std::unordered_set<face*>::iterator k = std::find(CH.begin(), CH.end(), *j);
      if(k != CH.end()){
        if(verbose) std::cout << "Deleting " << faces(*k) << " from Convex Hull\n";
        CH.erase(k);
      }
    }
    
    // Create faces and edges for edges in Horizon List
    std::vector<halfEdge*> horizon;
    std::vector<halfEdge*> delBuffer;
    std::vector<face*> delBufferF;
    obtainHorizon(v[i], facetList, horizon);
    halfEdge *before = new halfEdge;
    before->origin = v[i];
    v[i]->inc = before;
    for(int j = 0; j < horizon.size(); j++){
      halfEdge *e0 = new halfEdge, *e1 = new halfEdge, *e2;
      if(j == horizon.size()-1) e2 = v[i]->inc;
      else  e2 = new halfEdge;
      face* f_old = horizon[j]->twin->incident;
      horizon[j]->origin->inc = e1;
      if(j == horizon.size()-1) e0->origin = horizon[0]->origin;
      else  e0->origin = horizon[j+1]->origin;
      e1->origin = horizon[j]->origin;
      e2->origin = v[i];
      e0->twin = horizon[j]; horizon[j]->twin = e0;
      e1->twin = before; before->twin = e1;
      e0->next = e1; e1->next = e2; e2->next = e0;
      e0->prev = e2; e2->prev = e1; e1->prev = e0;
      std::cout << edges(horizon[0]);
      face *f = new face;
      e0->incident = e1->incident = e2->incident = f;
      f->boundary = e0;
      before = e2;

      if(verbose) std::cout << "Created new face " << faces(f) << " with " << coordinate(v[i]) << "\n";

      // If new face is coplanar with its adjacent face, merge the two
      face *g = e0->twin->incident;
      if(coplanar(f, g)){
        if(verbose) std::cout << "Merging new face with " << faces(g) << "\n";
        e0->twin->prev->next = e0->next;
        e0->twin->next->prev = e0->prev;
        e0->prev->next = e0->twin->next;
        e0->next->prev = e0->twin->prev;
        e1->incident = e2->incident = g;
        //Also delete e0 and e0->twin
        //delBuffer.push_back(e0);
        delete f;
        g->boundary = e1;
        
      } else {
      // Otherwise update conflict list and add the new face to Convex Hull
        if(verbose) std::cout << "Updating conflict graph for new face\n";
        // Get the list of vertices in conflict with the facet
        std::map<face*, std::unordered_set<vertex*> >::iterator it = faceConflict.find(f_old);
        delBufferF.push_back(f_old);
        if(it != faceConflict.end()){
          for(std::unordered_set<vertex*>::iterator cv = it->second.begin(); cv != it->second.end(); ++cv){
            // For each vertex originally in conflict with the old face, check for conflict with the new face
            if(*cv == v[i]) continue;
            if(isVisible(f, *cv)){
              if(verbose) std::cout << "Added " << faces(f) << " to CFD of " << coordinate(*cv) << "\n";
              faceConflict[f].insert(*cv);
              vertexConflict[*cv].insert(f);
            }
          }
        }
        it = faceConflict.find(g);
        if(it != faceConflict.end()){
          for(std::unordered_set<vertex*>::iterator cv = it->second.begin(); cv != it->second.end(); ++cv){
            if(*cv == v[i]) continue;
            if(isVisible(f, *cv)){
              faceConflict[f].insert(*cv);
              vertexConflict[*cv].insert(f);
            }
          }
        }
        CH.insert(f);
      }
      std::map<face*, std::unordered_set<vertex*> >::iterator it = faceConflict.find(f_old);
      for(std::unordered_set<vertex*>::iterator cv = it->second.begin(); cv != it->second.end(); ++cv){
        std::unordered_set<face*> &cvConflictList = vertexConflict[*cv];
        if(*cv == v[i]) continue;
        if(verbose) std::cout << "CFD for " << coordinate(*cv) << "\n";
        if(verbose) for(std::unordered_set<face*>::iterator ja = cvConflictList.begin(); ja != cvConflictList.end(); ++ja)
          std::cout << "\t" << faces(*ja) << "\n";
      }
    }
    deleteAll(delBufferF, v[i]);
    if(verbose) detail();
  }
}


void print(){
  std::set<vertex*> set;
  for(std::unordered_set<face*>::iterator it = CH.begin(); it != CH.end(); ++it){
    set.insert((*it)->boundary->origin);
    for(halfEdge* i = (*it)->boundary->next; i != (*it)->boundary; i = i->next)
      set.insert(i->origin);
  }
  std::cout << set.size() << std::endl;
  for(std::set<vertex*>::iterator it = set.begin(); it != set.end(); ++it)
    std::cout << (*it)->x << " " << (*it)->y << " " << (*it)->z << std::endl;
}

int main(){
  std::vector<vertex*> input;
  int n;
  
  std::cin >> n;
  for(int i = 0; i < n; i++){
    vertex *buff = new vertex;
    std::cin >> buff->x >> buff->y >> buff->z;
    buff->inc = nullptr;
    input.push_back(buff);
  }
 
  switch(makeTetrahedron(input, n)){
    case NONE:
    case LINEAR:
      //for(int i = 0; i < 2; i++) std::cout << input[i]->x << " " << input[i]->y << " " << input[i]->z << std::endl;
    case PLANAR:
      //PlanarHull(input);
      if(verbose) std::cout << "Tetrahedron does not exist\n";
      break;
    case TETRA:
      makeConflict(input, n);
      TetraHull(input, n);
      print();
      break;
 }
}
