#ifndef _MY_MESH_
#define _MY_MESH_


#include "Mesh\Vertex.h"
#include "Mesh\Edge.h"
#include "Mesh\Face.h"
#include "Mesh\HalfEdge.h"
#include "Mesh\BaseMesh.h"

#include "Mesh\boundary.h"
#include "Mesh\iterators.h"
#include "Parser\parser.h"
#include "Eigen\Dense"
#include <vector>
#ifndef M_PI
#define M_PI 3.141592653589793238
#endif

using Eigen::MatrixXd;
using namespace std;
namespace MeshLib
{
	class CMyVertex;
	class CMyEdge;
	class CMyFace;
	class CMyHalfEdge;

	
	class CMyVertex : public CVertex
	{
	public:
		CMyVertex() : m_rgb(1,1,1) {};
		~CMyVertex() {};

		void _from_string() ;
		void _to_string();

		CPoint & rgb() { return m_rgb; };
	protected:
		CPoint m_rgb;
	};

	inline void CMyVertex::_from_string()
	{
		CParser parser(m_string);
		for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
		{
			CToken * token = *iter;
			if (token->m_key == "uv") //CPoint2
			{
				token->m_value >> m_uv;
			}
			if (token->m_key == "rgb") // CPoint
			{
				token->m_value >> m_rgb;
			}
		}
	}

	inline void CMyVertex::_to_string()
	{
		CParser parser(m_string);
		parser._removeToken("uv");

		parser._toString(m_string);
		std::stringstream iss;

		iss << "uv=(" << m_uv[0] << " " << m_uv[1] << ")";

		if (m_string.length() > 0)
		{
			m_string += " ";
		}
		m_string += iss.str();
	}
	
	class CMyEdge : public CEdge
	{
	public:
		CMyEdge() :m_sharp(false) {};
		~CMyEdge() {};

		void _from_string();

		bool & sharp() { return m_sharp; };
	protected:
		bool m_sharp;
	};

	inline void CMyEdge::_from_string()
	{
		CParser parser(m_string);
		for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
		{
			CToken * token = *iter;
			if (token->m_key == "sharp") // bool
			{
				m_sharp = true;
			}
		}
	}

	class CMyFace : public CFace
	{
	public:

		CPoint & normal() { return m_normal; };
	protected:
		CPoint m_normal;
	};

	class CMyHalfEdge : public CHalfEdge
	{
	};

	class ValidPair
	{
	public:
		ValidPair(CVertex * v1, CVertex * v2, bool isedge = false)
		{
			m_vertex[0] = v1;
			m_vertex[1] = v2;
			is_edge = isedge;
		}
		~ValidPair();
		bool operator < (const ValidPair & vp)
		{
			return cost < vp.cost;
		}
		CVertex * & operator [] (int i)
		{
			assert(i > -1 && i < 2);
			return m_vertex[i];
		}
	private:
		CVertex * m_vertex[2];
		double cost;
		bool is_edge;
	};

	template<typename V, typename E, typename F, typename H>
	class MyMesh : public CBaseMesh<V, E, F, H>
	{
	public:
		typedef CBoundary<V, E, F, H>					CBoundary;
		typedef CLoop<V, E, F, H>						CLoop;

		typedef MeshVertexIterator<V, E, F, H>			MeshVertexIterator;
		typedef MeshEdgeIterator<V, E, F, H>			MeshEdgeIterator;
		typedef MeshFaceIterator<V, E, F, H>			MeshFaceIterator;
		typedef MeshHalfEdgeIterator<V, E, F, H>		MeshHalfEdgeIterator;

		typedef VertexVertexIterator<V, E, F, H>		VertexVertexIterator;
		typedef VertexEdgeIterator<V, E, F, H>			VertexEdgeIterator;
		typedef VertexFaceIterator<V, E, F, H>			VertexFaceIterator;
		typedef VertexInHalfedgeIterator<V, E, F, H>	VertexInHalfedgeIterator;
		typedef VertexOutHalfedgeIterator<V, E, F, H>	VertexOutHalfedgeIterator;

		typedef FaceVertexIterator<V, E, F, H>			FaceVertexIterator;
		typedef FaceEdgeIterator<V, E, F, H>			FaceEdgeIterator;
		typedef FaceHalfedgeIterator<V, E, F, H>		FaceHalfedgeIterator;

		void output_mesh_info();
		void test_iterator();
		vector<F *> faces_around_vertex(V * pV);
		MatrixXd plane_equation(F * pF);
		MatrixXd get_kp(MatrixXd m);
		MatrixXd compute_q(vector<F *> facesVec);
		MatrixXd compute_error(V * pV);
		double compute_cost(MatrixXd v, MatrixXd mQ);
		MatrixXd after_contraction(ValidPair vp);
	};

	typedef MyMesh<CMyVertex, CMyEdge, CMyFace, CMyHalfEdge> CMyMesh;

	template<typename V, typename E, typename F, typename H>
	void MeshLib::MyMesh<V, E, F, H>::output_mesh_info()
	{
		int nv = this->numVertices();
		int ne = this->numEdges();
		int nf = this->numFaces();

		std::cout << "#V=" << nv << "  ";
		std::cout << "#E=" << ne << "  ";
		std::cout << "#F=" << nf << "  ";

		int euler_char= nv - ne + nf;
		std::cout << "Euler's characteristic=" << euler_char << "  ";

		CBoundary boundary(this);
		std::vector<CLoop*> & loops = boundary.loops();
		int nb = loops.size();

		int genus = (2 - (euler_char + nb)) / 2;
		std::cout << "genus=" << genus << std::endl;
	}

	template<typename V, typename E, typename F, typename H>
	void MyMesh<V, E, F, H>::test_iterator()
	{
		for (MeshVertexIterator viter(this); !viter.end(); ++viter)
		{
			V * pV = *viter;
			// you can do something to the vertex here
			// ...

			for (VertexEdgeIterator veiter(pV); !veiter.end(); ++veiter)
			{
				E * pE = *veiter;
				// you can do something to the neighboring edges with CCW
				// ...
			}

			for (VertexFaceIterator vfiter(pV); !vfiter.end(); ++vfiter)
			{
				F * pF = *vfiter;
				// you can do something to the neighboring faces with CCW
				// ...
			}

			for (VertexInHalfedgeIterator vhiter(this, pV); !vhiter.end(); ++vhiter)
			{
				H * pH = *vhiter;
				// you can do something to the incoming halfedges with CCW
				// ...
			}
		}

		for (MeshEdgeIterator eiter(this); !eiter.end(); ++eiter)
		{
			E * pE = *eiter;
			// you can do something to the edge here
			// ...
		}

		for (MeshFaceIterator fiter(this); !fiter.end(); ++fiter)
		{
			F * pF = *fiter;
			// you can do something to the face here
			// ...
		}

		//there are some other iterators which you can find them in class MyMesh

		std::cout << "Iterators test OK.\n";
	}
	template<typename V, typename E, typename F, typename H>
	vector<F *> MyMesh<V, E, F, H>::faces_around_vertex(V * pV)
	{
		H * pH = pV->halfedge();
		H * pH_next = pH->next()->sym();
		vector<F *> facesVec;
		facesVec.push_back(pH->face());
		//traverse all faces, pV->halfedge() is a most ccw in-coming halfedge
		//thus we don't need traverse clwly
		while (pH_next != NULL && pH_next != pH)
		{
			facesVec.push_back(pH_next->face());
			pH_next = pH_next->next()->sym();
		}
		return facesVec;
	}

	template<typename V, typename E, typename F, typename H>
	MatrixXd MyMesh<V, E, F, H>::plane_equation(F * pF)
	{
		MatrixXd x(1, 4);
		MatrixXd A(3, 4);
		
		V * pV1 = pF->halfedge()->source();
		V * pV2 = pF->halfedge()->target();
		V * pV3 = pF->halfedge()->next()->target();
		//construct A matrix
		A << pV1->point()[0], pV1->point()[1], pV1->point()[2], 1,
			pV2->point()[0], pV2->point()[1], pV2->point()[2], 1,
			pV3->point()[0], pV3->point()[1], pV3->point()[2], 1;
		//solve linear equation Ax = B to get the plane equation
		x = A.lu().solve(MatrixXd::Zero(1, 4));
		//normalize x to let a^2 + b^2 + c^2 = 1
		//note that in eigen matrix starts at 0, not 1
		double sum = x(0, 0) * x(0, 0) + x(0, 1) * x(0, 1) + x(0, 2) * x(0, 2);
		x /= sqrt(sum);

		return x;
	}

	template<typename V, typename E, typename F, typename H>
	MatrixXd MyMesh<V, E, F, H>::get_kp(MatrixXd m)
	{
		return m * m.adjoint();
	}

	template<typename V, typename E, typename F, typename H>
	MatrixXd MyMesh<V, E, F, H>::compute_q(vector<F *> facesVec)
	{
		MatrixXd Q = MatrixXd::Zero(4, 4);
		for (vector::iterator it = facesVec.begin(); it != facesVec.end(); ++it)
		{
			Q += get_kp(plane_equation(*it));
		}
		return Q;
	}

	template<typename V, typename E, typename F, typename H>
	MatrixXd MyMesh<V, E, F, H>::compute_error(V * pV)
	{
		return compute_q(faces_around_vertex(pV));
	}
	
	/*
	cost = v' * Q * v
	*/
	template<typename V, typename E, typename F, typename H>
	double MyMesh<V, E, F, H>::compute_cost(MatrixXd v, MatrixXd mQ)
	{
		return v.adjoint() * mQ * v;
	}

	/*
	returns the position (x, y, z, 1) of a vertex after contraction
	param 1 is the Q matrix of v1 + v2
	*/
	template<typename V, typename E, typename F, typename H>
	MatrixXd MyMesh<V, E, F, H>::after_contraction(ValidPair vp)
	{
		MatrixXd mQ = compute_error(vp[0]) + compute_error(vp[1]);
		MatrixXd v(4, 1);
		MatrixXd A(4, 4);
		MatrixXd B(4, 1);
		//construct A
		A << mQ(0, 0), mQ(0, 1), mQ(0, 2), mQ(0, 3),
			mQ(1, 0), mQ(1, 1), mQ(1, 2), mQ(1, 3),
			mQ(2, 0), mQ(2, 1), mQ(2, 2), mQ(2, 3),
			0, 0, 0, 1;
		//if A is invertible
		if (A.determinant() != 0)
		{
			//construct B
			B << 0, 0, 0, 1;
			//solve the equation Ax = B
			v = A.lu().solve(B);
			return v
		}
		else //If the matrix is not invertible, then we can simply choose the optimal v
			//along the segment v1 v2.
		{
			double cost = 100000.0;
			float final_k = -1.0;
			for (float k = 0; k <= 1; k += 0.05)
			{
				v(0, 0) = (1 - k) * vp[0]->point()[0] + k * vp[1]->point()[0];
				v(1, 0) = (1 - k) * vp[0]->point()[1] + k * vp[1]->point()[1];
				v(2, 0) = (1 - k) * vp[0]->point()[2] + k * vp[1]->point()[2];
				v(3, 0) = 1;

				double tmp_cost = compute_cost(v, mQ);
				if (cost >= tmp_cost)
				{
					final_k = k;
					cost = tmp_cost;
				}
			}
			if (final_k != -1.0)
			{
				v(0, 0) = (1 - final_k) * vp[0]->point()[0] + final_k * vp[1]->point()[0];
				v(1, 0) = (1 - final_k) * vp[0]->point()[1] + final_k * vp[1]->point()[1];
				v(2, 0) = (1 - final_k) * vp[0]->point()[2] + final_k * vp[1]->point()[2];
				v(3, 0) = 1;
				return v;
			}
			else
			{
				cout << "error in k! \n";
				return NULL;
			}
		}
	}
}

#endif // !_MY_MESH_
