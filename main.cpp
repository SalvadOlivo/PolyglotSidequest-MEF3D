#include <iostream>
#include "math_tools.h"
#include "classes.h"
#include "tools.h"
#include "sel.h"

int main(int argc, char *argv[]){
    char filename[150];
    strcpy(filename,argv[1]);

    vector<Matrix> localKs;
    vector<Vector> localbs;
    Matrix K;
    Vector b;
    Vector T;

    cout << "IMPLEMENTACION DEL METODO DE LOS ELEMENTOS FINITOS\n"
         << "\t- TRANSFERENCIA DE CALOR\n" << "\t- 3 DIMENSIONES\n"
         << "\t- FUNCIONES DE FORMA LINEALES\n" << "\t- PESOS DE GALERKIN\n"
         << "\t- ELEMENTOS TETRAHEDROS\n"
         << "*********************************************************************************\n\n";

    mesh m;
    leerMallayCondiciones(m,filename);
    cout << "Datos obtenidos correctamente\n********************\n";

    crearSistemasLocales(m,localKs,localbs);
    //showKs(localKs); showbs(localbs);
    cout << "******************************\n";

    zeroes(K,m.getSize(NODES));
    zeroes(b,m.getSize(NODES));
    ensamblaje(m,localKs,localbs,K,b);
    //showMatrix(K); showVector(b);
    cout << "******************************\n";
    //cout << K.size() << " - "<<K.at(0).size()<<"\n";
    //cout << b.size() <<"\n";

    applyNeumann(m,b);
    //showMatrix(K); showVector(b);
    cout << "******************************\n";
    //cout << K.size() << " - "<<K.at(0).size()<<"\n";
    //cout << b.size() <<"\n";

    applyDirichlet(m,K,b);
    //showMatrix(K); showVector(b);
    cout << "******************************\n";
    //cout << K.size() << " - "<<K.at(0).size()<<"\n";
    //cout << b.size() <<"\n";

    zeroes(T,b.size());
    calculate(K,b,T);

    //cout << "La respuesta es: \n";
    //showVector(T);

    writeResults(m,T,filename);

    return 0;
}
