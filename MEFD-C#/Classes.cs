using System;
using System.Collections.Generic;
using System.Text;
using static MEF3D.Math_Tools;
using static MEF3D.Sel;
using static MEF3D.tools;

namespace MEF3D
{
    class Classes
    {
		public enum indicators
		{
			NOTHING
		}
		public enum lines
		{
			NOLINE,
			SINGLELINE,
			DOUBLELINE
		}
		public enum modes
		{
			NOMODE,
			INT_FLOAT,
			INT_FLOAT_FLOAT_FLOAT,
			INT_INT_INT_INT_INT
		}
		public enum parameters
		{
			THERMAL_CONDUCTIVITY,
			HEAT_SOURCE
		}
		public enum sizes
		{
			NODES,
			ELEMENTS,
			DIRICHLET,
			NEUMANN
		}

		public abstract class item
		{
			protected int id;
			protected float x;
			protected float y;
			protected float z;
			protected int node1;
			protected int node2;
			protected int node3;
			protected int node4;
			protected float value;
			public void setId(int identifier)
			{
				id = identifier;
			}

			public void setX(float x_coord)
			{
				x = x_coord;
			}

			public void setY(float y_coord)
			{
				y = y_coord;
			}

			public void setZ(float z_coord)
			{
				z = z_coord;
			}

			public void setNode1(int node_1)
			{
				node1 = node_1;
			}

			public void setNode2(int node_2)
			{
				node2 = node_2;
			}

			public void setNode3(int node_3)
			{
				node3 = node_3;
			}

			public void setNode4(int node_4)
			{
				node4 = node_4;
			}

			public void setValue(float value_to_assign)
			{
				value = value_to_assign;
			}

			public int getId()
			{
				return id;
			}

			public float getX()
			{
				return x;
			}

			public float getY()
			{
				return y;
			}

			public float getZ()
			{
				return z;
			}

			public int getNode1()
			{
				return node1;
			}

			public int getNode2()
			{
				return node2;
			}

			public int getNode3()
			{
				return node3;
			}

			public int getNode4()
			{
				return node4;
			}

			public float getValue()
			{
				return value;
			}

			public abstract void setValues(int a, float b, float c, float d, int e, int f, int g, int h, float i);

		}

		public class node : item
		{

			public void setValues(int a, float b, float c, float d, int e, int f, int g, int h, float i)
			{
				id = a;
				x = b;
				y = c;
				z = d;
			}

		}

		public class element : item
		{

			public void setValues(int a, float b, float c, float d, int e, int f, int g, int h, float i)
			{
				id = a;
				node1 = e;
				node2 = f;
				node3 = g;
				node4 = h;
			}

		}

		public class condition : item
		{


			public void setValues(int a, float b, float c, float d, int e, int f, int g, int h, float i)
			{
				node1 = e;
				value = i;
			}

		}

	public class mesh
		{
			private float[] parameters = new float[2];
			private int[] sizes = new int[4];
			private node[] node_list;
			private element[] element_list;
			private int[] indices_dirich;
			private condition[] dirichlet_list;
			private condition[] neumann_list;
			public void setParameters(float k, float Q)
			{
				parameters[THERMAL_CONDUCTIVITY] = k;
				parameters[HEAT_SOURCE] = Q;
			}
			public void setSizes(int nnodes, int neltos, int ndirich, int nneu)
			{
				sizes[NODES] = nnodes;
				sizes[ELEMENTS] = neltos;
				sizes[DIRICHLET] = ndirich;
				sizes[NEUMANN] = nneu;
			}
			public int getSize(int s)
			{
				return sizes[s];
			}
			public float getParameter(int p)
			{
				return parameters[p];
			}
			public void createData()
			{
				node_list = Arrays.InitializeWithDefaultInstances<node>(sizes[NODES]);
				element_list = Arrays.InitializeWithDefaultInstances<element>(sizes[ELEMENTS]);
				indices_dirich = new int[DIRICHLET];
				dirichlet_list = Arrays.InitializeWithDefaultInstances<condition>(sizes[DIRICHLET]);
				neumann_list = Arrays.InitializeWithDefaultInstances<condition>(sizes[NEUMANN]);
			}
			public node getNodes()
			{
				return node_list;
			}
			public element getElements()
			{
				return element_list;
			}
			//Acá hay que arreglar la logica con Dirichlet
			public int getDirichletIndices()
			{
				return indices_dirich;
			}
			public condition getDirichlet()
			{
				return dirichlet_list;
			}
			public condition getNeumann()
			{
				return neumann_list;
			}
			public node getNode(int i)
			{
				//esto da error
				return new node(node_list[i]);
			}
			public element getElement(int i)
			{
				//esto da error
				return new element(element_list[i]);
			}
			public condition getCondition(int i, int type)
			{
				if (type == DIRICHLET)
				{
					//esto da error
					return new condition(dirichlet_list[i]);
				}
				else
				{
					//todos estos dan error creo que por el constructor
					return new condition(neumann_list[i]);
				}
			}
		}




	}
}
