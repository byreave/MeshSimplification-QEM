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
#include <queue>
#include <algorithm>
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
		ValidPair(CVertex * v1, CVertex * v2, bool isedge = true)
		{
			m_vertex[0] = v1;
			m_vertex[1] = v2;
			is_edge = isedge;
		}
		ValidPair(){}
		~ValidPair(){}
		friend bool operator < (const ValidPair & vp1, const ValidPair & vp2)
		{
			return vp1.cost < vp2.cost;
		}
		friend bool operator > (const ValidPair & vp1, const ValidPair & vp2)
		{
			return vp1.cost > vp2.cost;
		}
		CVertex * & operator [] (int i)
		{
			assert(i > -1 && i < 2);
			return m_vertex[i];
		}
		void setCost(double co) { cost = co; }
		bool isEdge(){ return is_edge; }
		bool equal(ValidPair other)
		{
			return (m_vertex[0] == other.m_vertex[0] && m_vertex[1] == other.m_vertex[1]) || (m_vertex[0] == other.m_vertex[1] && m_vertex[1] == other.m_vertex[0]);
		}
		CVertex * & vertex(int i)
		{
			assert(i >= 0 && i < 2);
			return m_vertex[i];
		}
		void setVertex(int i, CVertex * v)
		{
			assert(i >= 0 && i < 2);
			m_vertex[i] = v;
		}
		void swap()
		{
			CVertex * tmp = m_vertex[0];
			m_vertex[0] = m_vertex[1];
			m_vertex[1] = tmp;
		}
		double getCost(){ return cost; }
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
		static vector<F *> faces_around_vertex(V * pV);
		static MatrixXd plane_equation(F * pF);
		static MatrixXd get_kp(MatrixXd m);
		static MatrixXd compute_q(vector<F *> facesVec);
		static MatrixXd compute_error(V * pV);
		static double compute_cost(MatrixXd v, MatrixXd mQ);
		static MatrixXd after_contraction(ValidPair vp);
		bool contract_valid_pair(ValidPair vp);
		void QEM(double threshold, int times);
	};

	typedef MyMesh<CMyVertex, CMyEdge, CMyFace, CMyHalfEdge> CMyMesh;

	class MyPriorityQueue{
	private:
		vector<ValidPair> data;

	public:
		void push(ValidPair vp){
			data.push_back(vp);
			push_heap(data.begin(), data.end(), greater<ValidPair>()); //min heap
		}

		void pop(){
			pop_heap(data.begin(), data.end(), greater<ValidPair>());
			data.pop_back();
		}

		ValidPair top() { return data.front(); }
		int size() { return data.size(); }
		bool empty() { return data.empty(); }
		bool hasPair(ValidPair * vp)
		{
			for (int i = 0; i < data.size(); ++i)
			{
				if (data[i].equal(*vp));
					return true;
			}
			return false;
		}
		void changeV1toV2(CVertex * v1, CVertex * v2)
		{
			//find all pair consists of v1
			for (int i = 0; i < data.size(); ++i)
			{
				if (data[i].vertex(0) == v1)
				{
					ValidPair * vp = new ValidPair(v2, data[i].vertex(1));
					for (vector<ValidPair>::iterator it = data.begin(); it != data.end();) //remove all duplicate paira
					{
						ValidPair itvp = *it;
						if (itvp.equal(*vp))
						{
							it = data.erase(it);
							make_heap(data.begin(), data.end(), greater<ValidPair>());
						}
						else
							it++;
					}
					delete vp;
					vp = NULL;
					data[i].setVertex(0, v2);
					MatrixXd v = CMyMesh::after_contraction(data[i]);
					data[i].setCost(CMyMesh::compute_cost(v, CMyMesh::compute_error((CMyVertex *)v2) + CMyMesh::compute_error((CMyVertex *)data[i].vertex(1))));
					make_heap(data.begin(), data.end(), greater<ValidPair>());
				}
				else if (data[i].vertex(1) == v1)
				{
					ValidPair * vp = new ValidPair(v2, data[i].vertex(0));
					for (vector<ValidPair>::iterator it = data.begin(); it != data.end();) //remove all duplicate paira
					{
						ValidPair itvp = *it;
						if (itvp.equal(*vp))
						{
							it = data.erase(it);
							make_heap(data.begin(), data.end(), greater<ValidPair>());
						}
						else
							it++;
					}
					delete vp;
					vp = NULL;
					data[i].setVertex(1, v2);
					MatrixXd v = CMyMesh::after_contraction(data[i]);
					data[i].setCost(CMyMesh::compute_cost(v, CMyMesh::compute_error((CMyVertex *)v2) + CMyMesh::compute_error((CMyVertex *)data[i].vertex(0))));
					make_heap(data.begin(), data.end(), greater<ValidPair>());
				}
				else
					continue;
			}
		}
	};

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
			//CHalfEdge * pH = edgeHalfedge(pE, 0);

		}

		for (MeshFaceIterator fiter(this); !fiter.end(); ++fiter)
		{
			F * pF = *fiter;
			cout << plane_equation(pF) << endl;
			break;
			// you can do something to the face here
			// ...
		}

		V * pV = idVertex(264);
		std::cout << pV->edges().size();

		H * pH = (H *)pV->halfedge();
		H * pH_next = (H *)pH->he_next()->he_sym();
		cout << pH->source()->id() << endl;
		while (pH_next != NULL && pH_next != pH)
		{
			cout << pH_next->source()->id() << endl;
			pH_next = (H *)pH_next->he_next()->he_sym();
		}
		//there are some other iterators which you can find them in class MyMesh

		
		std::cout << "Iterators test OK.\n";
	}
	template<typename V, typename E, typename F, typename H>
	vector<F *> MyMesh<V, E, F, H>::faces_around_vertex(V * pV)
	{
		H * pH = (H *)pV->halfedge();
		H * pH_next;
		
		vector<F *> facesVec;
		if (pH != NULL)
			pH_next = (H *)pH->he_next()->he_sym();
		else
			return facesVec;
		facesVec.push_back((F *)pH->face());
		//traverse all faces, pV->halfedge() is a most ccw in-coming halfedge
		//thus we don't need traverse clwly
		while (pH_next != NULL && pH_next != pH)
		{
			facesVec.push_back((F *)pH_next->face());
			pH_next = (H *)pH_next->he_next()->he_sym();
		}
		return facesVec;
	}

	template<typename V, typename E, typename F, typename H>
	MatrixXd MyMesh<V, E, F, H>::plane_equation(F * pF)
	{
		MatrixXd x(4, 1);
		
		V * pV1 = (V *)pF->halfedge()->source();
		V * pV2 = (V *)pF->halfedge()->target();
		V * pV3 = (V *)pF->halfedge()->he_next()->target();

		double a = (pV2->point()[1] - pV1->point()[1]) * (pV3->point()[2] - pV1->point()[2]) - (pV3->point()[1] - pV1->point()[1]) * (pV2->point()[2] - pV1->point()[2]);
		double b = (pV2->point()[2] - pV1->point()[2]) * (pV3->point()[0] - pV1->point()[0]) - (pV3->point()[2] - pV1->point()[2]) * (pV2->point()[0] - pV1->point()[0]);
		double c = (pV2->point()[0] - pV1->point()[0]) * (pV3->point()[1] - pV1->point()[1]) - (pV3->point()[0] - pV1->point()[0]) * (pV2->point()[1] - pV1->point()[1]);
		double d = -a * pV1->point()[0] - b * pV1->point()[1] - c * pV1->point()[2];
		
		x << a, b, c, d;
		//normalize x to let a^2 + b^2 + c^2 = 1
		//note that in eigen matrix starts at 0, not 1
		double sum = x(0, 0) * x(0, 0) + x(1, 0) * x(1, 0) + x(1, 0) * x(1, 0);
		if (sum == 0)
			return MatrixXd::Zero(4, 1);
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
		for (vector<F *>::iterator it = facesVec.begin(); it != facesVec.end(); ++it)
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
		return (v.adjoint() * mQ * v).determinant();
	}

	/*
	returns the position (x, y, z, 1) of a vertex after contraction
	*/
	template<typename V, typename E, typename F, typename H>
	MatrixXd MyMesh<V, E, F, H>::after_contraction(ValidPair vp)
	{
		MatrixXd Q1 = compute_error((V *)vp[0]);
		MatrixXd Q2 = compute_error((V *)vp[1]);

		MatrixXd mQ = Q1 + Q2;
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
			cout << "SOLVE" << endl;
			//construct B
			B << 0, 0, 0, 1;
			//solve the equation Ax = B
			v = A.lu().solve(B);
			//cout << "v=" << v << endl;
			if (_isnan(v(0, 0)) != 0)
			{
				cout << "Eception occured" << vp[0]->id() << ' ' << vp[1]->id() << endl;
				cout << Q1 << endl;
				
				assert(1 < 0);
			}
			return v;
		}
		else //If the matrix is not invertible, then we can simply choose the optimal v
			//along the segment v1 v2.
		{
			cout << "INTER" << endl;

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
				return v;
			}
		}
	}
	template<typename V, typename E, typename F, typename H>
	bool MyMesh<V, E, F, H>::contract_valid_pair(ValidPair vp)
	{
		MatrixXd v = after_contraction(vp); //get position after contraction
		vp[0]->point()[0] = v(0, 0);
		vp[0]->point()[1] = v(1, 0);
		vp[0]->point()[2] = v(2, 0);

		

		cout <<"after con :"<< v << endl;
		if (vp.isEdge())
		{
			///new///
			if (vp[0]->boundary() && !vp[1]->boundary())
			{
				E * e = (E *)vertexEdge((V *)vp[0], (V *)vp[1]);

				//find the halfedge with source being vp[0] and target being vp[1]
				H * h;
				if (e->halfedge(0)->source() == vp[0])
					h = (H *)e->halfedge(0);
				else if (e->halfedge(1) != NULL)
					h = (H *)e->halfedge(1);

				//vertex
				//m_edges
				if (vp[0]->id() < vp[1]->id())
					vp[0]->edges().remove(e);
				H * pH1 = vp[1]->halfedge();
				if (pH1->source()->id() > vp[0]->id())
					vp[0]->edges().push_back(pH1->edge());
				H * pH1_next = pH1->he_next()->he_sym();
				while (pH1_next != pH1)
				{
					if (pH1_next->source()->id() > vp[0]->id())
						vp[0]->edges().push_back(pH1_next->edge());
					pH1_next = pH1_next->he_next()->he_sym();
				}
				//halfedge
				if (vp[0]->halfedge() == h->prev())
					vp[0]->halfedge() = h->he_next()->he_sym();
				if (h->->he_sym()->he_next()->target()->halfedge() == h->he_sym()->he_next())
					h->he_next()->target()->halfedge() = h->he_sym()->he_prev()->he_sym();
				//target
				pH1_next = pH1->he_next()->he_sym();
				pH1->target() = vp[0];
				while (pH1_next != pH1)
				{
					pH1_next->target() = vp[0];
				}

				m_verts.remove((V*)vp[1]);
				std::map<int, V*>::iterator viter = m_map_vert.find(vp[1]->id());
				if (viter != m_map_vert.end())
				{
					m_map_vert.erase(viter);
				}
				

				//edge
				//halfedges
				if (h->he_prev()->edge()->halfedge(0) == h->he_prev())
					h->he_prev()->edge()->halfedge(0) = h->he_next()->he_sym();
				else
					h->he_prev()->edge()->halfedge(1) = h->he_next()->he_sym();
				h->he_next()->he_sym()->edge() = h->he_prev()->edge();
				if (h->he_sym()->he_next()->edge()->halfedge(0) == h->he_sym()->he_next())
					h->he_sym()->he_next()->edge()->halfedge(0) = h->he_sym()->he_prev()->he_sym();
				else
					h->he_sym()->he_next()->edge()->halfedge(1) = h->he_sym()->he_prev()->he_sym();
				h->he_sym()->he_prev()->he_sym()->edge() = h->he_sym()->he_next()->edge();
			}



			///new///
			E * e = (E *)vertexEdge((V *)vp[0], (V *)vp[1]);
			
			//find the halfedge with source being vp[0] and target being vp[1]
			H * h;
			if (e->halfedge(0)->source() == vp[0])
				h = (H *)e->halfedge(0);
			else if (e->halfedge(1) != NULL)
				h = (H *)e->halfedge(1);
			else //former operation makes sure there is a halfedge vp[0]->vp[1]
			{
				cout << "error in halfedge\n";
				return false;
			}
			if (vp[1]->boundary())
			{
				if (!vp[0]->boundary())
					vp[0]->boundary() = true;
				else //in case a single edge comes up after contraction
					return false;
			}
			H * h_prev = (H *)h->he_prev();
			H * h_next = (H *)h->he_next();

			//remove vp[1]
			//m_edges in vp[1] to vp[0]
			if (vp[0]->id() < vp[1]->id())
				vp[0]->edges().remove(e);

			if (h_next->target()->id() < vp[1]->id())
				h_next->target()->edges().remove(h_next->edge());

			if (h->he_sym() != NULL)
			{
				if (h->he_sym()->he_next()->target()->id() < vp[1]->id())
					h->he_sym()->he_next()->target()->edges().remove(h->he_sym()->he_prev()->edge());

			}
			H * hvp1 = (H *)vp[1]->halfedge();
			H * hvp1_next = (H *)hvp1->he_next()->he_sym();
			if (hvp1->source() != vp[0] && hvp1->source() != h_next->target())
			{
				if (h->he_sym() == NULL && hvp1->source()->id() > vp[0]->id()) //boundary
					vp[0]->edges().push_back(hvp1->edge());
				else if (h->he_sym() != NULL && hvp1->source() != h->he_sym()->he_next()->target() && hvp1->source()->id() > vp[0]->id())
					vp[0]->edges().push_back(hvp1->edge());
				
			}
			while (hvp1_next != NULL && hvp1_next != hvp1)
			{
				if (hvp1_next->source() != vp[0] && hvp1_next->source() != h_next->target())
				{
					if (h->he_sym() == NULL && hvp1_next->source()->id() > vp[0]->id()) //boundary
						vp[0]->edges().push_back(hvp1_next->edge());
					else if (h->he_sym() != NULL && hvp1_next->source() != h->he_sym()->he_next()->target() && hvp1_next->source()->id() > vp[0]->id())
						vp[0]->edges().push_back(hvp1_next->edge());
				}
				hvp1_next = (H *)hvp1_next->he_next()->he_sym();
			}

			//set all in and out halfedges of vp[1] to vp[0]
			H * h_in_vp1 = (H *)vp[1]->halfedge();
			H * h_out_vp1 = (H *)h_in_vp1->he_sym();
			H * h_out_vp1_next = h_out_vp1 == NULL ? NULL : (H *)h_out_vp1->he_next()->he_sym();
			H * h_in_vp1_next = (H *)h_in_vp1->he_next()->he_sym();
			h_in_vp1->target() = vp[0];
			while (h_in_vp1_next != h_in_vp1 && h_in_vp1_next != NULL)
			{
				h_in_vp1_next->target() = vp[0];
				h_in_vp1_next = (H *)h_in_vp1_next->he_next()->he_sym();
			}
			if (h_out_vp1 != NULL)
				h_out_vp1->source() = vp[0];
			while (h_out_vp1_next != h_out_vp1 && h_out_vp1_next != NULL)
			{
				h_out_vp1_next->source() = vp[0];
				h_out_vp1_next = (H *)h_out_vp1_next->he_next()->he_sym();
			}
			std::map<int, V*>::iterator viter = m_map_vert.find(vp[1]->id());
			if (viter != m_map_vert.end())
			{
				m_map_vert.erase(viter);
			}
			m_verts.remove((V *)vp[1]);//delete vp[1]; don't delete the memory because we need maintain the heap, see more in QEM()

			//edges
			H * h_prev_sym = (H *)h_prev->he_sym();
			if (h_prev->edge()->halfedge(0) == h_prev && h_prev->edge()->halfedge(1) != NULL)
			{
				h_prev->edge()->halfedge(0) = h_prev->edge()->halfedge(1);
				h_prev->edge()->halfedge(1) = h_next->he_sym();

			}
			else if (h_prev->edge()->halfedge(0) != h_prev)
				h_prev->edge()->halfedge(1) = h_next->he_sym();
			else if (h_prev->edge()->halfedge(0) == h_prev && h_prev->edge()->halfedge(1) == NULL)
				h_prev->edge()->halfedge(0) = h_next->he_sym();
			if (h_next->he_sym() != NULL) //not boundary
				h_next->he_sym()->edge() = h_prev->edge();
			m_edges.remove((E *)h_next->edge());

			if (h_next->target()->halfedge() == h_next)
				h_next->target()->halfedge() = h_prev_sym;
			if (h_prev->target()->halfedge() == h_prev && h_next->he_sym() != NULL)
				h_prev->target()->halfedge() = h_next->he_sym();
			else if (h_prev->target()->halfedge() == h_prev && h_next->he_sym() == NULL)
			{
				h_prev->target()->halfedge() =h_prev_sym->he_prev();
			}
			else if (h->he_sym() != NULL && h_prev->target()->halfedge() == h->he_sym() && h_next->he_sym() != NULL)
			{
				h_prev->target()->halfedge() = h_next->he_sym();
			}
			else if (h->he_sym() != NULL && h_prev->target()->halfedge() == h->he_sym() && h_next->he_sym() == NULL)
			{
				h_prev->target()->halfedge() = h_prev_sym->he_prev();
			}
			
			delete h_next->edge();
			h_next->edge() = NULL;
			delete h_prev;
			delete h_next;
			h_next = NULL;
			h_prev = NULL;
			

			if (h->he_sym() != NULL) // edge for contraction isn't on the boundary
			{
				h = (H *)h->he_sym();
				h_prev = (H *)h->he_prev();
				h_next = (H *)h->he_next();
				H * h_next_sym = (H *)h_next->he_sym();
				//edges same as above
				if (h_next->edge()->halfedge(0) == h_next && h_next->edge()->halfedge(1) != NULL)
				{
					h_next->edge()->halfedge(0) = h_next->edge()->halfedge(1);
					h_next->edge()->halfedge(1) = h_prev->he_sym();
				}
				else if (h_next->edge()->halfedge(0) != h_next)
					h_next->edge()->halfedge(1) = h_prev->he_sym();
				else if (h_next->edge()->halfedge(0) == h_next && h_next->edge()->halfedge(1) == NULL)
					h_next->edge()->halfedge(0) = h_prev->he_sym();
				if (h_prev->he_sym() != NULL)
					h_prev->he_sym()->edge() = h_next->edge();
				m_edges.remove((E *)h_prev->edge());

				if (h_next->target()->halfedge() == h_next && h_prev->he_sym() != NULL)
					h_next->target()->halfedge() = h_prev->he_sym();
				else if (h_next->target()->halfedge() == h_next && h_prev->he_sym() == NULL)
					h_next->target()->halfedge() = h_next_sym->he_prev();
				delete h_prev->edge();
				h_prev->edge() = NULL;
				delete h_next;
				delete h_prev;
				h_prev = NULL;
				h_next = NULL;
			}

			//delete faces beside e
			std::map<int, F *>::iterator fiter = m_map_face.find(h->face()->id());
			if (fiter != m_map_face.end())
			{
				m_map_face.erase(fiter);
				m_faces.remove((F *)h->face());

				delete h->face();
			}
			
			if (h->he_sym() != NULL)
			{
				h = (H *)h->he_sym();
				fiter = m_map_face.find(h->face()->id());
				if (fiter != m_map_face.end())
				{
					m_map_face.erase(fiter);
					m_faces.remove((F *)h->face());

					delete h->face();
				}
				
				delete h->he_sym();
			}
			//remove edge between vp[0] and vp[1]
			m_edges.remove((E *)h->edge());
			delete h->edge();
			h->edge() = NULL;
			delete h;
		}
		else //if vp[0] vp[1] don't share an edge vp[1] must be on boundary
		{
			H * hvp1 = (H *)(vp[1]->halfedge());//most ccwly in-coming halfedge must be on boundary
			H * hvp1_next = (H *)(hvp1->he_next()->he_sym());
			//m_edges
			if (hvp1->source()->id() < vp[0]->id())
				vp[0]->edges().push_back((E *)hvp1->edge());
			//halfedges
			hvp1->target() = vp[0];
			while (hvp1_next != NULL)
			{
				if (hvp1_next->source()->id() < vp[0]->id())
					vp[0]->edges().push_back((E *)hvp1_next->edge());
				hvp1_next->target() = vp[0];
				if (hvp1_next->he_sym() != NULL)
				{
					hvp1_next->he_sym()->source() = vp[0];
				}
				hvp1_next = (H *)hvp1_next->he_next()->he_sym();
			}
			std::map<int, V*>::iterator viter = m_map_vert.find(vp[1]->id());
			if (viter != m_map_vert.end())
			{
				m_map_vert.erase(viter);
			}
			m_verts.remove((V *)vp[1]);//delete vp[1];
		}
		return true;
	}
	template<typename V, typename E, typename F, typename H>
	void MyMesh<V, E, F, H>::QEM(double threshold, int times)
	{
		MyPriorityQueue heap;
		for (MeshEdgeIterator eiter(this); !eiter.end(); ++eiter)//all valid pairs on edge
		{
			E * pE = *eiter;
			V * v1 = (V *)pE->halfedge(0)->source();
			V * v2 = (V *)pE->halfedge(0)->target();
			ValidPair * vp = new ValidPair(v1, v2);
			MatrixXd v = after_contraction(*vp);
			vp->setCost(compute_cost(v, compute_error(v1) + compute_error(v2)));
				
			heap.push(*vp);
		}

		//add valid pairs that are not on edge
		for (MeshVertexIterator viter1(this); !viter1.end(); ++viter1)
		{
			V * pV = *viter1;
			for (MeshVertexIterator viter2(this); !viter2.end(); ++viter2)
			{
				V * pV2 = *viter2;
				if (vertexEdge(pV, pV2) == NULL)
				{
					ValidPair *vp = new ValidPair(pV, pV2);
					if (!heap.hasPair(vp))
					{
						MatrixXd mQ = compute_error(pV) + compute_error(pV2);
						double cost = compute_cost(after_contraction(*vp), mQ);
						if (cost <= threshold)
						{
							vp->setCost(cost);
							heap.push(*vp);
						}
					}
				}
			}
		}
		ValidPair vp;
		//contraction of valid pairs with minimum cost for (times) times
		for (int i = 0; i < times && !heap.empty(); ++i)
		{
			vp = heap.top();
			if (contract_valid_pair(vp))
			{
				//cout << vp.getCost() << endl;
				cout << vp[0]->id() << ' ' << vp[1]->id() << endl;
				heap.pop();
				//change all changed pairs of the deleted vertex to the pair of remaining vertex and delete duplicate pairs
				heap.changeV1toV2(vp[1], vp[0]);
				//delete vp[1];
			}
			else
			{
				heap.pop();
			}
		}
	}
}

#endif // !_MY_MESH_
