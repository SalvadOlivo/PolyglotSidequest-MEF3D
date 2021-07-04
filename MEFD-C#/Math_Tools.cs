using System;
using System.Collections.Generic;
using System.Text;
using static MEF3D.Classes;
using static MEF3D.Sel;
using static MEF3D.tools;

namespace MEF3D
{
    class Math_Tools
    {
		private void zeroes(List<List<float>> M, int n)
		{
			for (int i = 0; i < n; i++)
			{
				List<float> row = new List<float>(n);
				M.Add(row);
			}
		}

		private void zeroes(List<List<float>> M, int n, int m)
		{
			for (int i = 0; i < n; i++)
			{
				List<float> row = new List<float>(m);
				M.Add(row);
			}
		}

		private void zeroes(List<float> v, int n)
		{
			for (int i = 0; i < n; i++)
			{
				v.Add(0.0F);
			}
		}

		private void copyMatrix(List<List<float>> A, List<List<float>> copy)
		{
			zeroes(copy, A.Count);
			for (int i = 0; i < A.Count; i++)
			{
				for (int j = 0; j < A[0].Count; j++)
				{
					copy[i][j] = A[i][j];
				}
			}
		}

		private float calculateMember(int i, int j, int r, List<List<float>> A, List<List<float>> B)
		{
			float member = 0F;
			for (int k = 0; k < r; k++)
			{
				member += A[i][k] * B[k][j];
			}
			return member;
		}

		private List<List<float>> productMatrixMatrix(List<List<float>> A, List<List<float>> B, int n, int r, int m)
		{
			List<List<float>> R = new List<List<float>>();

			zeroes(R, n, m);
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					R[i][j] = calculateMember(i, j, r, new List<List<float>>(A), new List<List<float>>(B));
				}
			}

			return new List<List<float>>(R);
		}

		private void productMatrixVector(List<List<float>> A, List<float> v, List<float> R)
		{
			for (int f = 0; f < A.Count; f++)
			{
				float cell = 0.0F;
				for (int c = 0; c < v.Count; c++)
				{
					cell += A[f][c] * v[c];
				}
				R[f] += cell;
			}
		}

		private void productRealMatrix(float real, List<List<float>> M, List<List<float>> R)
		{
			zeroes(R, M.Count);
			for (int i = 0; i < M.Count; i++)
			{
				for (int j = 0; j < M[0].Count; j++)
				{
					R[i][j] = real * M[i][j];
				}
			}
		}

		private void getMinor(List<List<float>> M, int i, int j)
		{
			//cout << "Calculando menor ("<<i+1<<","<<j+1<<")...\n";
			M.RemoveAt(i);
			for (int i = 0; i < M.Count; i++)
			{
				//no recuerdo cual es el equivalente para el metodo "erase" 
				M[i].erase(M[i].GetEnumerator() + j);
			}
		}

		private float determinant(List<List<float>> M)
		{
			if (M.Count == 1)
			{
				return M[0][0];
			}
			else
			{
				float det = 0.0F;
				for (int i = 0; i < M[0].Count; i++)
				{
					List<List<float>> minor = new List<List<float>>();
					copyMatrix(new List<List<float>>(M), minor);
					getMinor(minor, 0, i);
                    (void)(det += Math.Pow(-1, i) * M[0][i] * determinant(new List<List<float>>(minor)));
				}
				return det;
			}
		}

		private void cofactors(List<List<float>> M, List<List<float>> Cof)
		{
			zeroes(Cof, M.Count);
			for (int i = 0; i < M.Count; i++)
			{
				for (int j = 0; j < M[0].Count; j++)
				{
					//cout << "Calculando cofactor ("<<i+1<<","<<j+1<<")...\n";
					List<List<float>> minor = new List<List<float>>();
					copyMatrix(new List<List<float>>(M), minor);
					getMinor(minor, i, j);
					Cof[i][j] = (float)(Math.Pow(-1, i + j) * determinant(new List<List<float>>(minor)));
				}
			}
		}

		private void transpose(List<List<float>> M, List<List<float>> T)
		{
			zeroes(T, M[0].Count, M.Count);
			for (int i = 0; i < M.Count; i++)
			{
				for (int j = 0; j < M[0].Count; j++)
				{
					T[j][i] = M[i][j];
				}
			}
		}

		private void inverseMatrix(List<List<float>> M, List<List<float>> Minv)
		{
			Console.Write("Iniciando calculo de inversa...\n");
			List<List<float>> Cof = new List<List<float>>();
			List<List<float>> Adj = new List<List<float>>();
			Console.Write("Calculo de determinante...\n");
			float det = determinant(new List<List<float>>(M));
			if (det == 0F)
			{
				Environment.Exit(1);
			}
			Console.Write("Iniciando calculo de cofactores...\n");
			cofactors(new List<List<float>>(M), Cof);
			Console.Write("Calculo de adjunta...\n");
			transpose(new List<List<float>>(Cof), Adj);
			Console.Write("Calculo de inversa...\n");
			productRealMatrix(1 / det, new List<List<float>>(Adj), Minv);
		}
	}
}
