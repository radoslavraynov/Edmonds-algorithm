// Radoslav Zlatev, rzlatev@uni-bonn.de
// Programming Assignment, 26.11.2008

#ifndef COMBINATORIALOPTIMIZATION_RADOSLAVZLATEV_CLASSES2
#define COMBINATORIALOPTIMIZATION_RADOSLAVZLATEV_CLASSES2

#include <vector>
#include <iostream>
#include <fstream>

class Vertex;
class Graph;

/*
 * the class vertex implements (for convenience) the needed functions
 * on a vertex -- the set of neighbours, the \mu, \phi, and \ro vertices
 * and the scanned function
 */
class Vertex {
private:
	unsigned int _id;
	std::vector<Vertex*> _nbhrs;
	
	Vertex *_mu, *_fi, *_ro;
	bool _scanned;

public:
/* constructor taking only an integer index of the vertex */
	Vertex(unsigned int id) {
		_id = id;
		_mu = _fi = _ro = this;
		_scanned = false;
	};
/* empty destructor */
	~Vertex(void) { };
	
/* getters */
	Vertex* mu(void) { return _mu; };
	Vertex* fi(void) { return _fi; };
	Vertex* ro(void) { return _ro; };
	Vertex* getNbhr(unsigned int i) { if (i<_nbhrs.size()) return _nbhrs[i]; else return NULL; };
	unsigned int numNbhrs(void) { return _nbhrs.size(); };
	unsigned int id(void) { return _id; };

/* setters */
	void setMu(Vertex *v) { _mu = v; };
	void setFi(Vertex *v) { _fi = v; };
	void setRo(Vertex *v) { _ro = v; };
	void setNbhr(Vertex *v) { _nbhrs.push_back(v); };
	void setScanned(bool boolval) { _scanned = boolval; };
	
/* property methods */
	bool scanned(void) { return _scanned; };
	bool outer(void) { return (this->mu() == this) || (this->mu()->fi() != this->mu()); };
	bool inner(void) { return (this->mu()->fi() == this->mu()) && (this->fi() != this); };
	bool outOfForest(void) { return !(this->inner() || this->outer()); };
};

/*
 * the Graph class implements the graph itself as container of the vertices
 * it has an empty constructor and is initialized by the step1() procedure
 * the steps stepN() as methods represent the respective Nth step of the
 * algorithm; each step N, thus, either call step N+1 with the assumed parameters
 * or calles a previous step. At termination, the matching is printed as required
 */
class Graph {
private:
	std::vector<Vertex*> G;
	unsigned int n;
	
public:
	Graph(void) { };
	~Graph(void) {
		for (unsigned int i=0; i<n; i++)
			delete G[i];
	};
	
	std::vector<Vertex*> makePath(unsigned int r, Vertex* u = NULL) {
		std::vector<Vertex*> path;
		Vertex *v, *w = G[r];
		bool valid = true;
		
		// this case should never occur in the excecution of the algorithm
		if (w == u)
			return std::vector<Vertex*>(1,w);
		
		for (unsigned int i=0; valid; i++) {
			v = w;
			path.push_back(v);
			if (i%2)
				w = v->fi();
			else
				w = v->mu();
			
			// check if path will remain a path after adding w
			for (unsigned int j=0; j<=i; j++)
				if ((path[j]==w) || (u==w)) {
					valid = false;
					break;
				}
		}
		
		if (u==w)
			path.push_back(w);
		
		return path;
	};
	
	/*
	 * read the graph off from the specification file and puts a pointer
	 * to an instance of Vertex for every vertex from the graph; puts
	 * the neighbours in the respective container. NOTE: works only if
	 * there is no empty line at the end of the file
	 */
	void step1(char* filename) {
		int vi, vj;
	
		std::ifstream gspec(filename);
	
		gspec >> n;
	
		for (unsigned int i=0; i<n; i++) {
			Vertex* v = new Vertex(i);
			G.push_back(v);
		}
	
		while (!gspec.eof()) {
			gspec >> vi >> vj;
			G[vi]->setNbhr(G[vj]);
			G[vj]->setNbhr(G[vi]);
		}
	};
	
	/* as in the algorithm given in class */
	void step2(void) {
		int x = -1;
		for (unsigned int i=0; i<n; i++)
			if (G[i]->outer() && !(G[i]->scanned()) && x==-1) {
				x = i;
				break;
			}
	
		if (x == -1)
			this->printMatching();
		else
			this->step3(x);
	};
	
	/* as in the algorithm given in class */
	void step3(unsigned int x) {
		int y = -1;
		for (unsigned int i=0; i<G[x]->numNbhrs(); i++)
			if ((G[x]->getNbhr(i)->outOfForest()) || (G[x]->getNbhr(i)->outer() && (G[x]->getNbhr(i)->ro() != G[x]->ro()))) {
				// such a vertex y was found, so we break and go to step 4.
				y = G[x]->getNbhr(i)->id();
				break;
			}

		if (y == -1) {
			// there was no such y
			G[x]->setScanned(true);
			this->step2();
		}
		else
			this->step4(x,y);
	};
	
	/* as in the algorithm given in class */
	void step4(unsigned int x, unsigned int y) {
		if (G[y]->outOfForest()) {
			G[y]->setFi(G[x]);
			this->step3(x);
		}
		else
 			this->step5(x,y);
	};
	
	/* as in the algorithm given in class */
	void step5(unsigned int x, unsigned int y) {
		std::vector<Vertex*> Px(this->makePath(x)), Py(this->makePath(y));
		
		/*
		 * Note that for the algorithm to work two paths are considered vertex disjoint
		 * if they do not *any* vertices and not if they do not common inner vertices.
		 */
		for (unsigned int i=0; i<Px.size(); i++)
			for (unsigned int j=0; j<Py.size(); j++)
				if (Px[i] == Py[j]) {
					/*
					 * Here note that we do not need to consider the case
					 * that the paths are not disjont but for no vertex ro(v)=v,
					 * because the condition that y is outer assures that if they
					 * meet then they meet also in some root of blossom
					 */
					if (Px[i]->ro() == Px[i]) {
						this->step6(x,y,Px[i]->id());
						return;
					}
				}

		for (unsigned int i=1; i<Px.size(); i+=2) {
			Px[i]->fi()->setMu(Px[i]);
			Px[i]->setMu(Px[i]->fi());
		}
		for (unsigned int i=1; i<Py.size(); i+=2) {
			Py[i]->fi()->setMu(Py[i]);
			Py[i]->setMu(Py[i]->fi());
		}
		G[x]->setMu(G[y]);
		G[y]->setMu(G[x]);
		
		for (unsigned int j=0; j<n; j++) {
			G[j]->setFi(G[j]);
			G[j]->setRo(G[j]);
			G[j]->setScanned(false);
		}
		
		this->step2();
	};
	
	/* as in the algorithm given in class */
	void step6(unsigned int x, unsigned int y, unsigned int r) {
		std::vector<Vertex*> Px(this->makePath(x,G[r])), Py(this->makePath(y,G[r]));
		
		for (unsigned int i=1; i<Px.size(); i+=2)
			if ((Px[i]->fi()->ro() != Px[i]) && (Px[i]->ro() != Px[i]))
				Px[i]->fi()->setFi(Px[i]);
		
		for (unsigned int i=1; i<Py.size(); i+=2)
			if ((Py[i]->fi()->ro() != Py[i]) && (Py[i]->ro() != Py[i]))
				Py[i]->fi()->setFi(Py[i]);
		
		if (G[x]->ro() != G[r])
			G[x]->setFi(G[y]);
		if (G[y]->ro() != G[r])
			G[y]->setFi(G[x]);
		
		for (unsigned int j=0; j<n; j++) {
			bool tocheck = true;
			for (unsigned int k=0; k<Px.size(); k++)
				if (G[j]->ro() == Px[k]) {
					G[j]->setRo(G[r]);
					tocheck = false;
					break;
				}
			for (unsigned int k=0; tocheck && (k<Py.size()); k++)
				if (G[j]->ro() == Py[k]) {
					G[j]->setRo(G[r]);
					break;
				}
		}
		
		this->step3(x);
	};
	
	/* print the matching in the form: one edge "i j" with indeces i < j on every line */
	void printMatching() {
		int j=0;
		for (unsigned int i=0; i<n; i++)
			if (i < G[i]->mu()->id()) {
				std::cout << i << " " << G[i]->mu()->id() << std::endl;
				j++;
			}
		std::cout << std::endl << j << std::endl;
	};
};

#endif // COMBINATORIALOPTIMIZATION_RADOSLAVZLATEV_CLASSES2
