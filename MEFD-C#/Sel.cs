using System;
using System.Collections.Generic;
using System.Text;
using static MEF3D.Classes;
using static MEF3D.Math_Tools;
using static MEF3D.tools;
using System.Drawing.Drawing2D;
using System.Numerics;


//BUSCAR LIBRERIA PARA MANEJAR LAS MATRICES Y VECTORES
namespace MEF3D 
{
    class Sel
    {
		private void showMatrix(Matrix K)
		{
			for (int i = 0; i < K.at(0).size(); i++)
			{
				Console.Write("[\t");
				for (int j = 0; j < K.size(); j++)
				{
					Console.Write(K.at(i).at(j));
					Console.Write("\t");
				}
				Console.Write("]\n");
			}
		}

		private void showKs(List<Matrix> Ks)
		{
			for (int i = 0; i < Ks.Count; i++)
			{
				Console.Write("K del elemento ");
				Console.Write(i + 1);
				Console.Write(":\n");
				showMatrix(new List<Matrix>(Ks[i]));
				Console.Write("*************************************\n");
			}
		}

		private void showVector(Vector b)
		{
			Console.Write("[\t");
			for (int i = 0; i < b.size(); i++)
			{
				Console.Write(b.at(i));
				Console.Write("\t");
			}
			Console.Write("]\n");
		}

		private void showbs(List<Vector> bs)
		{
			for (int i = 0; i < bs.Count; i++)
			{
				Console.Write("b del elemento ");
				Console.Write(i + 1);
				Console.Write(":\n");
				showVector(new List<Vector>(bs[i]));
				Console.Write("*************************************\n");
			}
		}

		private float calculateLocalD(int ind, mesh m)
		{
			float D;
			float a;
			float b;
			float c;
			float d;
			float e;
			float f;
			float g;
			float h;
			float i;

			element el = m.getElement(ind);

			node n1 = m.getNode(el.getNode1() - 1);
			node n2 = m.getNode(el.getNode2() - 1);
			node n3 = m.getNode(el.getNode3() - 1);
			node n4 = m.getNode(el.getNode4() - 1);

			a = n2.getX() - n1.getX();
			b = n2.getY() - n1.getY();
			c = n2.getZ() - n1.getZ();
			d = n3.getX() - n1.getX();
			e = n3.getY() - n1.getY();
			f = n3.getZ() - n1.getZ();
			g = n4.getX() - n1.getX();
			h = n4.getY() - n1.getY();
			i = n4.getZ() - n1.getZ();
			//Se calcula el determinante de esta matriz utilizando
			//la Regla de Sarrus.
			D = a * e * i + d * h * c + g * b * f - g * e * c - a * h * f - d * b * i;

			return D;
		}

		private float calculateLocalVolume(int ind, mesh m)
		{
			//Se utiliza la siguiente fÃ³rmula:
			//      Dados los 4 puntos vÃ©rtices del tetrahedro A, B, C, D.
			//      Nos anclamos en A y calculamos los 3 vectores:
			//              V1 = B - A
			//              V2 = C - A
			//              V3 = D - A
			//      Luego el volumen es:
			//              V = (1/6)*det(  [ V1' ; V2' ; V3' ]  )

			float V;
			float a;
			float b;
			float c;
			float d;
			float e;
			float f;
			float g;
			float h;
			float i;
			element el = m.getElement(ind);
			node n1 = m.getNode(el.getNode1() - 1);
			node n2 = m.getNode(el.getNode2() - 1);
			node n3 = m.getNode(el.getNode3() - 1);
			node n4 = m.getNode(el.getNode4() - 1);

			a = n2.getX() - n1.getX();
			b = n2.getY() - n1.getY();
			c = n2.getZ() - n1.getZ();
			d = n3.getX() - n1.getX();
			e = n3.getY() - n1.getY();
			f = n3.getZ() - n1.getZ();
			g = n4.getX() - n1.getX();
			h = n4.getY() - n1.getY();
			i = n4.getZ() - n1.getZ();
			//Para el determinante se usa la Regla de Sarrus.
			V = (1.0 / 6.0) * (a * e * i + d * h * c + g * b * f - g * e * c - a * h * f - d * b * i);

			return V;
		}

		private float ab_ij(float ai, float aj, float a1, float bi, float bj, float b1)
		{
			return (ai - a1) * (bj - b1) - (aj - a1) * (bi - b1);
		}

		private void calculateLocalA(int i, Matrix A, mesh m)
		{
			element e = m.getElement(i);
			node n1 = m.getNode(e.getNode1() - 1);
			node n2 = m.getNode(e.getNode2() - 1);
			node n3 = m.getNode(e.getNode3() - 1);
			node n4 = m.getNode(e.getNode4() - 1);

			A.at(0).at(0) = ab_ij(n3.getY(), n4.getY(), n1.getY(), n3.getZ(), n4.getZ(), n1.getZ());
			A.at(0).at(1) = ab_ij(n4.getY(), n2.getY(), n1.getY(), n4.getZ(), n2.getZ(), n1.getZ());
			A.at(0).at(2) = ab_ij(n2.getY(), n3.getY(), n1.getY(), n2.getZ(), n3.getZ(), n1.getZ());
			A.at(1).at(0) = ab_ij(n4.getX(), n3.getX(), n1.getX(), n4.getZ(), n3.getZ(), n1.getZ());
			A.at(1).at(1) = ab_ij(n2.getX(), n4.getX(), n1.getX(), n2.getZ(), n4.getZ(), n1.getZ());
			A.at(1).at(2) = ab_ij(n3.getX(), n2.getX(), n1.getX(), n3.getZ(), n2.getZ(), n1.getZ());
			A.at(2).at(0) = ab_ij(n3.getX(), n4.getX(), n1.getX(), n3.getY(), n4.getY(), n1.getY());
			A.at(2).at(1) = ab_ij(n4.getX(), n2.getX(), n1.getX(), n4.getY(), n2.getY(), n1.getY());
			A.at(2).at(2) = ab_ij(n2.getX(), n3.getX(), n1.getX(), n2.getY(), n3.getY(), n1.getY());
		}

		private void calculateB(Matrix B)
		{
			B.at(0).at(0) = -1;
			B.at(0).at(1) = 1;
			B.at(0).at(2) = 0;
			B.at(0).at(3) = 0;
			B.at(1).at(0) = -1;
			B.at(1).at(1) = 0;
			B.at(1).at(2) = 1;
			B.at(1).at(3) = 0;
			B.at(2).at(0) = -1;
			B.at(2).at(1) = 0;
			B.at(2).at(2) = 0;
			B.at(2).at(3) = 1;
		}

		private Matrix createLocalK(int element, mesh m)
		{
			// K = (k*Ve/D^2)Bt*At*A*B := K_4x4
			float D;
			float Ve;
			float k = m.getParameter(THERMAL_CONDUCTIVITY);
			Matrix K = new Matrix();
			Matrix A = new Matrix();
			Matrix B = new Matrix();
			Matrix Bt = new Matrix();
			Matrix At = new Matrix();

			D = calculateLocalD(element, m);
			Ve = calculateLocalVolume(element, m);

			zeroes(A, 3);
			zeroes(B, 3, 4);
			calculateLocalA(element, A, m);
			calculateB(B);
			transpose(A, At);
			transpose(B, Bt);

			productRealMatrix(k * Ve / (D * D), productMatrixMatrix(Bt, productMatrixMatrix(At, productMatrixMatrix(A, B, 3, 3, 4), 3, 3, 4), 4, 3, 4), K);

			return new Matrix(K);
		}

		private float calculateLocalJ(int ind, mesh m)
		{
			float J;
			float a;
			float b;
			float c;
			float d;
			float e;
			float f;
			float g;
			float h;
			float i;

			element el = m.getElement(ind);

			node n1 = m.getNode(el.getNode1() - 1);
			node n2 = m.getNode(el.getNode2() - 1);
			node n3 = m.getNode(el.getNode3() - 1);
			node n4 = m.getNode(el.getNode4() - 1);

			a = n2.getX() - n1.getX();
			b = n3.getX() - n1.getX();
			c = n4.getX() - n1.getX();
			d = n2.getY() - n1.getY();
			e = n3.getY() - n1.getY();
			f = n4.getY() - n1.getY();
			g = n2.getZ() - n1.getZ();
			h = n3.getZ() - n1.getZ();
			i = n4.getZ() - n1.getZ();
			//Se calcula el determinante de esta matriz utilizando
			//la Regla de Sarrus.
			J = a * e * i + d * h * c + g * b * f - g * e * c - a * h * f - d * b * i;

			return J;
		}

		private Vector createLocalb(int element, mesh m)
		{
			Vector b = new Vector();

			float Q = m.getParameter(HEAT_SOURCE);
			float J;
			float b_i;
			J = calculateLocalJ(element, new mesh(m));

			b_i = Q * J / 24.0;
			b.push_back(b_i);
			b.push_back(b_i);
			b.push_back(b_i);
			b.push_back(b_i);

			return new Vector(b);
		}

		private void crearSistemasLocales(mesh m, List<Matrix> localKs, List<Vector> localbs)
		{
			for (int i = 0; i < m.getSize(ELEMENTS); i++)
			{
				localKs.Add(createLocalK(i, m));
				localbs.Add(createLocalb(i, m));
			}
		}

		private void assemblyK(element e, Matrix localK, Matrix K)
		{
			int index1 = e.getNode1() - 1;
			int index2 = e.getNode2() - 1;
			int index3 = e.getNode3() - 1;
			int index4 = e.getNode4() - 1;

			K.at(index1).at(index1) += localK.at(0).at(0);
			K.at(index1).at(index2) += localK.at(0).at(1);
			K.at(index1).at(index3) += localK.at(0).at(2);
			K.at(index1).at(index4) += localK.at(0).at(3);
			K.at(index2).at(index1) += localK.at(1).at(0);
			K.at(index2).at(index2) += localK.at(1).at(1);
			K.at(index2).at(index3) += localK.at(1).at(2);
			K.at(index2).at(index4) += localK.at(1).at(3);
			K.at(index3).at(index1) += localK.at(2).at(0);
			K.at(index3).at(index2) += localK.at(2).at(1);
			K.at(index3).at(index3) += localK.at(2).at(2);
			K.at(index3).at(index4) += localK.at(2).at(3);
			K.at(index4).at(index1) += localK.at(3).at(0);
			K.at(index4).at(index2) += localK.at(3).at(1);
			K.at(index4).at(index3) += localK.at(3).at(2);
			K.at(index4).at(index4) += localK.at(3).at(3);
		}
		private void assemblyb(element e, Vector localb, Vector b)
		{
			int index1 = e.getNode1() - 1;
			int index2 = e.getNode2() - 1;
			int index3 = e.getNode3() - 1;
			int index4 = e.getNode4() - 1;

			b.at(index1) += localb.at(0);
			b.at(index2) += localb.at(1);
			b.at(index3) += localb.at(2);
			b.at(index4) += localb.at(3);
		}

		private void ensamblaje(mesh m, List<Matrix> localKs, List<Vector> localbs, Matrix K, Vector b)
		{
			for (int i = 0; i < m.getSize(ELEMENTS); i++)
			{
				element e = m.getElement(i);
				assemblyK(e, localKs[i], K);
				assemblyb(new element(e), new List<Vector>(localbs[i]), b);
			}
		}

		private void applyNeumann(mesh m, Vector b)
		{
			for (int i = 0; i < m.getSize(NEUMANN); i++)
			{
				condition c = m.getCondition(i, NEUMANN);
				b.at(c.getNode1() - 1) += c.getValue();
			}
		}

		private void applyDirichlet(mesh m, Matrix K, Vector b)
		{
			for (int i = 0; i < m.getSize(DIRICHLET); i++)
			{
				condition c = m.getCondition(i, DIRICHLET);
				int index = c.getNode1() - 1;

				K.erase(K.begin() + index);
				b.erase(b.begin() + index);

				for (int row = 0; row < K.size(); row++)
				{
					float cell = K.at(row).at(index);
					K.at(row).erase(K.at(row).begin() + index);
					b.at(row) += -1 * c.getValue() * cell;
				}
			}
		}

		private void calculate(Matrix K, Vector b, Vector T)
		{
			Console.Write("Iniciando calculo de respuesta...\n");
			Matrix Kinv = new Matrix();
			Console.Write("Calculo de inversa...\n");
			inverseMatrix(K, Kinv);
			Console.Write("Calculo de respuesta...\n");
			productMatrixVector(Kinv, b, T);
		}


	}
}
