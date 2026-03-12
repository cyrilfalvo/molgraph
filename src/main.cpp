
#include "constants.h"
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <fstream>
#include "nauty.h"
#include "gtools.h"

using namespace std;

void gen_topology(int natoms, int* types, double* coords, int* topology) {
  
  for(int i=0;i<natoms*natoms;i++) { topology[i] = 0; } 
  
  for(int i=0;i<natoms;i++) {
    for(int j=i+1;j<natoms;j++) {
      double dx = coords[i*3+0]-coords[j*3+0];
      double dy = coords[i*3+1]-coords[j*3+1];
      double dz = coords[i*3+2]-coords[j*3+2];
      double dr2 = dx*dx+dy*dy+dz*dz;
   
      double rmax2;

      if((types[i] == 6)&&(types[j] == 6)) { rmax2 = rcc2; }
      if((types[i] == 1)&&(types[j] == 6)) { rmax2 = rch2; }
      if((types[i] == 6)&&(types[j] == 1)) { rmax2 = rch2; }
      if((types[i] == 1)&&(types[j] == 1)) { rmax2 = rhh2; }

      if(dr2<=rmax2) { topology[i*natoms+j] = 1; topology[j*natoms+i] = 1; }
      //printf("%d %d %d\n",i,j,topology[i*natoms+j]);
    }
  }
  
}

void fixgroups(int i, int* groups, int natoms, int* topology) {
  for(int j=0;j<natoms;j++) {
    if(topology[i*natoms+j]>0) {
      if(groups[j] == 0) {
        groups[j] = groups[i];
	fixgroups(j,groups,natoms,topology);
      } else if (groups[j] != groups[i] ) {
            printf("ERROR HERE\n"); 
            printf("%d %d %d %d %d\n",i,j,topology[i*natoms+j],groups[i],groups[j]);
            exit(0);
      } 
    }
  }
}


int gen_nfrags(int natoms, int* topology) {
  int nfrags = 0;
  int* groups = new int[natoms];
  for(int i=0;i<natoms;i++) { groups[i] = 0; }
  for(int i=0;i<natoms;i++) {
    if(groups[i] == 0) {
      nfrags++;
      groups[i] = nfrags;
      fixgroups(i,groups,natoms,topology);
    }
  }
  delete [] groups;
  return nfrags;
}

int compute_n3(int natoms,int* type, int* topology) {
  int n3=0;
  for(int iatom=0;iatom<natoms;iatom++) {
    if(type[iatom]==1) {
      int nneighbors = 0;
      for(int jatom=0;jatom<natoms;jatom++) {
        if(type[jatom]==1) {
	  nneighbors += topology[iatom*natoms+jatom];
        }
      }   
      if(nneighbors==3) {
        n3 += 1;
      }
    }
  }
  return n3;
}


int main(int argc, char* argv[]) {

  string xyzfile;

  switch(argc) {
    case 2 : { 
      xyzfile = argv[1];
      break;
    }
    default : { 
      printf("usage: ./molgraph input.xyz\n");
      return EXIT_FAILURE;
    }
  }
  

  // Open the input file
  
  FILE *fp;
  fp = fopen(xyzfile.c_str(), "r");
  if (fp == NULL) {
      printf("can't open the input file\n");
      return EXIT_FAILURE;
  }

  int natoms;
  double *coords;
  int *types;
  int* topology;
  int *indextype;
  int ntypes;

  // read the number of atoms
  fscanf(fp, "%d", &natoms);
  fgetc(fp);
  
  // allocate arrays
  types = (int*)malloc(natoms*sizeof(int));
  coords = (double*)malloc(3*natoms*sizeof(double));
  topology = (int*)malloc(natoms*natoms*sizeof(int));

  // skip a line
  char buffer[256];
  fgets(buffer, sizeof(buffer), fp);

  // read the data
  char atype[3];
  for (int i=0;i<natoms;i++) {
    fscanf(fp, "%3s %lf %lf %lf",
             atype,&coords[3*i],&coords[3*i+1],&coords[3*i+2]);
    if(strcmp(atype,"H")==0) { types[i] = 1;}
    else if(strcmp(atype,"C")==0) { types[i] = 6; }
    else {
      printf("unknown atom type %s\n",atype);
      return EXIT_FAILURE;
    } 
  }

  // close the file
  fclose(fp);
  
  // sort the atoms
  for(int i=0;i<natoms;i++) {
    for(int j=i+1;j<natoms;j++) {
      if(types[i] < types[j] ) { 
        int itmp = types[i];
        types[i] = types[j];
        types[j] = itmp;
        for(int k=0;k<3;j++) {
          double dtmp = coords[3*i+k]; 
          coords[3*i+k] = coords[3*j+k];
	  coords[3*j+k] = dtmp;
	}
      }
    }
  }
  
  // count the number of types of atoms  
  ntypes = 1;
  for(int i=1;i<natoms;i++) {
    if(types[i] != types[i-1]) {
      ntypes++;
    }
  }

  printf("# found %d types of atoms\n",ntypes);

  indextype = (int*)malloc((ntypes+1)*sizeof(int));

  indextype[0] = 0;
  int itype=0;
  for(int i=1;i<natoms;i++) {
    if(types[i] != types[i-1]) {
      itype++;
      indextype[itype] = i;
    }
  }
  indextype[ntypes] = natoms;

  printf("# molecule: ");
  for(int i=0;i<ntypes;i++) {
    string atype;
    if(types[indextype[i]] == 1) {
      atype = "H";
    }	else if (types[indextype[i]] == 6) {
      atype = "C";
    } else {
        printf("unknown atom type\n");
        return EXIT_FAILURE;
    }
    printf("%s%d",atype.c_str(),indextype[i+1]-indextype[i]);
  }
  printf("\n");    
    
  gen_topology(natoms,types,coords,topology);
  int nfrags = gen_nfrags(natoms, topology);

  printf("# found %d fragment(s)\n",nfrags);


  DYNALLSTAT(int,lab,lab_sz);
  DYNALLSTAT(int,ptn,ptn_sz);
  DYNALLSTAT(int,orbits,orbits_sz);
  DYNALLSTAT(int,map,map_sz);
  DYNALLSTAT(graph,g,g_sz);
  DYNALLSTAT(graph,cg,cg_sz);
  static DEFAULTOPTIONS_GRAPH(options);
  statsblk stats;
  int n,m;
  //size_t k;
  //Select option for canonical labelling 
  options.getcanon = TRUE;
  options.defaultptn = FALSE;

  n = natoms;
  m = SETWORDSNEEDED(n);
  nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
  
  DYNALLOC1(int,lab,lab_sz,n,"malloc");
  DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
  DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
  DYNALLOC1(int,map,map_sz,n,"malloc");
  DYNALLOC2(graph,g,g_sz,n,m,"malloc");
  DYNALLOC2(graph,cg,cg_sz,n,m,"malloc");


  char* cg_str;

  for(int i=0;i<natoms;i++) { lab[i] = i; }
  for(int j=0;j<ntypes;j++) {
    for(int i=indextype[j];i<indextype[j+1]-1;i++) {
      ptn[i] = 1;
    }
    ptn[indextype[j+1]-1] = 0;
  }
  
  //for(int i=0;i<natoms;i++) {
  //  printf("%d %d %d\n",i,lab[i],ptn[i]);
  //}

  //Now make the graph
  EMPTYGRAPH(g,m,n);
  for(int i=0;i<n;i++) {
    for(int j=i+1;j<n;j++) {
      if(topology[i*n+j] == 1) {
        ADDONEEDGE(g,i,j,m);
      }  
    }
  }

  densenauty(g,lab,ptn,orbits,&options,&stats,m,n,cg);
  cg_str = (char*)ntos6(cg,m,n);
  printf("%s",cg_str);  
 
  DYNFREE(lab,lab_sz);
  DYNFREE(ptn,ptn_sz);
  DYNFREE(orbits,orbits_sz);
  DYNFREE(map,map_sz);
  DYNFREE(g,g_sz);
  DYNFREE(cg,cg_sz);

  
  free(indextype);
  free(types);
  free(coords);  
  free(topology);

  return EXIT_SUCCESS;

}
