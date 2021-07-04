#include <fstream>
#include "string.h"

void obtenerDatos(istream &file,int nlines,int n,int mode,item* item_list){
    string line;
    file >> line;
    if(nlines==DOUBLELINE) file >> line;

    for(int i=0;i<n;i++){
        switch(mode){
        case INT_FLOAT:
            int e0; float r0;
            file >> e0 >> r0;
            item_list[i].setValues(NOTHING,NOTHING,NOTHING,NOTHING,e0,NOTHING,NOTHING,NOTHING,r0);
            break;
        case INT_FLOAT_FLOAT_FLOAT:
            int e; float r,rr,rrr;
            file >> e >> r >> rr >> rrr;
            item_list[i].setValues(e,r,rr,rrr,NOTHING,NOTHING,NOTHING,NOTHING,NOTHING);
            break;
        case INT_INT_INT_INT_INT:
            int e1,e2,e3,e4,e5;
            file >> e1 >> e2 >> e3 >> e4 >> e5;
            item_list[i].setValues(e1,NOTHING,NOTHING,NOTHING,e2,e3,e4,e5,NOTHING);
            break;
        }
    }
}

void correctConditions(int n,condition *list,int *indices){
    for(int i=0;i<n;i++)
        indices[i] = list[i].getNode1();

    for(int i=0;i<n-1;i++){
        int pivot = list[i].getNode1();
        for(int j=i;j<n;j++)
            //Si la condición actual corresponde a un nodo posterior al nodo eliminado por
            //aplicar la condición anterior, se debe actualizar su posición.
            if(list[j].getNode1()>pivot)
                list[j].setNode1(list[j].getNode1()-1);
    }
}

void addExtension(char *newfilename,char *filename,char *extension){
    int ori_length = strlen(filename);
    int ext_length = strlen(extension);
    int i;
    for(i=0;i<ori_length;i++)
        newfilename[i] = filename[i];
    for(i=0;i<ext_length;i++)
        newfilename[ori_length+i] = extension[i];
    newfilename[ori_length+i] = '\0';
}

void leerMallayCondiciones(mesh &m,char *filename){
    char inputfilename[150];
    ifstream file;
    float k,Q;
    int nnodes,neltos,ndirich,nneu;

    addExtension(inputfilename,filename,".dat");
    file.open(inputfilename);

    file >> k >> Q;
    //cout << "k y Q: "<<k<<" y "<<Q<<"\n";
    file >> nnodes >> neltos >> ndirich >> nneu;
    //cout << "sizes: "<<nnodes<<" y "<<neltos<<" y "<<ndirich<<" y "<<nneu<<"\n";

    m.setParameters(k,Q);
    m.setSizes(nnodes,neltos,ndirich,nneu);
    m.createData();

    obtenerDatos(file,SINGLELINE,nnodes,INT_FLOAT_FLOAT_FLOAT,m.getNodes());
    obtenerDatos(file,DOUBLELINE,neltos,INT_INT_INT_INT_INT,m.getElements());
    obtenerDatos(file,DOUBLELINE,ndirich,INT_FLOAT,m.getDirichlet());
    obtenerDatos(file,DOUBLELINE,nneu,INT_FLOAT,m.getNeumann());

    file.close();

    //Se corrigen los índices en base a las filas que serán eliminadas
    //luego de aplicar las condiciones de Dirichlet
    correctConditions(ndirich,m.getDirichlet(),m.getDirichletIndices());
}

bool findIndex(int v, int s, int *arr){
    for(int i=0;i<s;i++)
        if(arr[i]==v) return true;
    return false;
}

void writeResults(mesh m,Vector T,char *filename){
    char outputfilename[150];
    int *dirich_indices = m.getDirichletIndices();
    condition *dirich = m.getDirichlet();
    ofstream file;

    addExtension(outputfilename,filename,".post.res");
    file.open(outputfilename);

    file << "GiD Post Results File 1.0\n";
    file << "Result \"Temperature\" \"Load Case 1\" 1 Scalar OnNodes\nComponentNames \"T\"\nValues\n";

    int Tpos = 0;
    int Dpos = 0;
    int n = m.getSize(NODES);
    int nd = m.getSize(DIRICHLET);
    for(int i=0;i<n;i++){
        if(findIndex(i+1,nd,dirich_indices)){
            file << i+1 << " " << dirich[Dpos].getValue() << "\n";
            Dpos++;
        }else{
            file << i+1 << " " << T.at(Tpos) << "\n";
            Tpos++;
        }
    }

    file << "End values\n";

    file.close();
}
