using System;
using System.Collections.Generic;
using static MEF3D.Classes;
using static MEF3D.Math_Tools;
using static MEF3D.Sel;
using static MEF3D.tools;
using System.Drawing.Drawing2D;
using System.Numerics;


//acá falta arreglar casi todo


namespace MEF3D
{
    class MainP
    {
        static void Main(string[] args)
        {
			string filename = new string(new char[150]);
			filename = args[1];

			List<Matrix> localKs = new List<Matrix>();
			List<Vector> localbs = new List<Vector>();
			Matrix K = new Matrix();
			Vector b = new Vector();
			Vector T = new Vector();

			Console.Write("IMPLEMENTACION DEL METODO DE LOS ELEMENTOS FINITOS\n");
			Console.Write("\t- MOVIMIENTO DE VIGA\n");
			Console.Write("\t- 3 DIMENSIONES\n");
			Console.Write("\t- FUNCIONES DE FORMA LINEALES\n");
			Console.Write("\t- PESOS DE GALERKIN\n");
			Console.Write("\t- ELEMENTOS TETRAHEDROS\n");
			Console.Write("*********************************************************************************\n\n");

			mesh m = new mesh();
			leerMallayCondiciones(m, filename);
			Console.Write("Datos obtenidos correctamente\n********************\n");

			crearSistemasLocales(m, localKs, localbs);
			//showKs(localKs); showbs(localbs);
			Console.Write("******************************\n");

			zeroes(K, m.getSize(NODES));
			zeroes(b, m.getSize(NODES));
			ensamblaje(m, localKs, localbs, K, b);
			//showMatrix(K); showVector(b);
			Console.Write("******************************\n");
			//cout << K.size() << " - "<<K.at(0).size()<<"\n";
			//cout << b.size() <<"\n";

			applyNeumann(m, b);
			//showMatrix(K); showVector(b);
			Console.Write("******************************\n");
			//cout << K.size() << " - "<<K.at(0).size()<<"\n";
			//cout << b.size() <<"\n";

			applyDirichlet(m, K, b);
			//showMatrix(K); showVector(b);
			Console.Write("******************************\n");
			//cout << K.size() << " - "<<K.at(0).size()<<"\n";
			//cout << b.size() <<"\n";

			zeroes(T, b.size());
			calculate(K, b, T);

			//cout << "La respuesta es: \n";
			//showVector(T);

			writeResults(m, T, filename);


		}
	}
}
