#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <set>
#include <random>
#include <chrono>
#include <sstream>
#include <SFML/Graphics.hpp>

using namespace std;

const double EPS = 1e-9;

struct Point {
    double x, y;
    Point(double x = 0, double y = 0) : x(x), y(y) {}

    bool operator<(const Point& p) const {
        if (fabs(x - p.x) > EPS) return x < p.x;
        return y < p.y - EPS;
    }

    bool operator==(const Point& p) const {
        return fabs(x - p.x) < EPS && fabs(y - p.y) < EPS;
    }

    Point operator-(const Point& p) const {
        return Point(x - p.x, y - p.y);
    }

    Point operator+(const Point& p) const {
        return Point(x + p.x, y + p.y);
    }

    Point operator*(double s) const {
        return Point(x * s, y * s);
    }

    double cross(const Point& p) const {
        return x * p.y - y * p.x;
    }

    double dot(const Point& p) const {
        return x * p.x + y * p.y;
    }

    double len2() const {
        return x * x + y * y;
    }

    double dist2(const Point& p) const {
        return (x - p.x)*(x - p.x) + (y - p.y)*(y - p.y);
    }
};

struct Circle {
    Point center;
    double radius;
    Circle() : center(0, 0), radius(0) {}
    Circle(Point c, double r) : center(c), radius(r) {}
    Circle(const Point& p1, const Point& p2) {
        center = Point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2);
        radius = sqrt(center.dist2(p1));
    }
    Circle(const Point& p1, const Point& p2, const Point& p3) {
        double d = 2 * (p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y));
        if (fabs(d) < EPS) {
            *this = Circle(p1, p2);
            return;
        }
        double ux = ((p1.x*p1.x + p1.y*p1.y) * (p2.y - p3.y) +
                     (p2.x*p2.x + p2.y*p2.y) * (p3.y - p1.y) +
                     (p3.x*p3.x + p3.y*p3.y) * (p1.y - p2.y)) / d;
        double uy = -((p1.x*p1.x + p1.y*p1.y) * (p2.x - p3.x) +
                      (p2.x*p2.x + p2.y*p2.y) * (p3.x - p1.x) +
                      (p3.x*p3.x + p3.y*p3.y) * (p1.x - p2.x)) / d;
        center = Point(ux, uy);
        radius = sqrt(center.dist2(p1));
    }

    bool contains(const Point& p) const {
        return center.dist2(p) <= radius * radius + EPS;
    }
};

struct Vertex;
struct Face;

struct HalfEdge {
    int id;
    Vertex* origin;
    HalfEdge* twin;
    HalfEdge* next;
    HalfEdge* prev;
    Face* face;

    HalfEdge() : id(-1), origin(nullptr), twin(nullptr),
                 next(nullptr), prev(nullptr), face(nullptr) {}
};

struct Vertex {
    int id;
    Point point;
    HalfEdge* incident;
    Vertex() : id(-1), incident(nullptr) {}
    Vertex(int id, Point p) : id(id), point(p), incident(nullptr) {}
};

struct Face {
    int id;
    HalfEdge* outer;
    vector<HalfEdge*> inner;

    Face() : id(-1), outer(nullptr) {}
};

struct DCEL {
    vector<Vertex*> vertices;
    vector<HalfEdge*> halfEdges;
    vector<Face*> faces;
    int nextVertexId;
    int nextHalfEdgeId;
    int nextFaceId;

    DCEL() : nextVertexId(0), nextHalfEdgeId(0), nextFaceId(0) {}

    ~DCEL() {
        for (size_t i = 0; i < vertices.size(); i++) delete vertices[i];
        for (size_t i = 0; i < halfEdges.size(); i++) delete halfEdges[i];
        for (size_t i = 0; i < faces.size(); i++) delete faces[i];
    }

    Vertex* addVertex(Point p) {
        Vertex* v = new Vertex(nextVertexId++, p);
        vertices.push_back(v);
        return v;
    }

    HalfEdge* addHalfEdge() {
        HalfEdge* he = new HalfEdge();
        he->id = nextHalfEdgeId++;
        halfEdges.push_back(he);
        return he;
    }

    Face* addFace() {
        Face* f = new Face();
        f->id = nextFaceId++;
        faces.push_back(f);
        return f;
    }

    void addEdge(Vertex* v1, Vertex* v2) {
        if (!v1 || !v2) return;

        HalfEdge* he1 = addHalfEdge();
        HalfEdge* he2 = addHalfEdge();

        he1->origin = v1;
        he2->origin = v2;
        he1->twin = he2;
        he2->twin = he1;

        he1->next = he2;
        he2->prev = he1;
        he2->next = he1;
        he1->prev = he2;

        if (!v1->incident) v1->incident = he1;
        if (!v2->incident) v2->incident = he2;
    }
};

struct Triangle {
    Point a, b, c;
    int ia, ib, ic;
    Circle circumcircle;

    Triangle() : ia(-1), ib(-1), ic(-1) {}

    Triangle(Point a, Point b, Point c, int ia, int ib, int ic)
        : a(a), b(b), c(c), ia(ia), ib(ib), ic(ic) {
        updateCircumcircle();
    }

    void updateCircumcircle() {
        circumcircle = Circle(a, b, c);
    }

    bool containsInCircumcircle(const Point& p) const {
        return circumcircle.contains(p);
    }

    bool hasVertex(int idx) const {
        return ia == idx || ib == idx || ic == idx;
    }

    bool hasEdge(int idx1, int idx2) const {
        if (idx1 < 0 || idx2 < 0) return false;
        return (ia == idx1 && ib == idx2) || (ia == idx2 && ib == idx1) ||
               (ib == idx1 && ic == idx2) || (ib == idx2 && ic == idx1) ||
               (ic == idx1 && ia == idx2) || (ic == idx2 && ia == idx1);
    }

    bool hasEdgePoints(const Point& p1, const Point& p2) const {
        return (a == p1 && b == p2) || (a == p2 && b == p1) ||
               (b == p1 && c == p2) || (b == p2 && c == p1) ||
               (c == p1 && a == p2) || (c == p2 && a == p1);
    }

    int getThirdVertex(int idx1, int idx2) const {
        if (ia != idx1 && ia != idx2) return ia;
        if (ib != idx1 && ib != idx2) return ib;
        if (ic != idx1 && ic != idx2) return ic;
        return -1;
    }
};

class MinimumEnclosingCircle {
private:
    vector<Point> points;
    vector<int> indices;

public:
    MinimumEnclosingCircle(const vector<Point>& pts) : points(pts) {
        indices.resize(points.size());
        for (size_t i = 0; i < points.size(); i++) indices[i] = i;
    }

    Circle compute() {
        if (points.empty()) return Circle(Point(0, 0), 0);
        if (points.size() == 1) return Circle(points[0], 0);

        shuffle(indices.begin(), indices.end(),
                mt19937(chrono::steady_clock::now().time_since_epoch().count()));

        Circle c(points[indices[0]], 0);

        for (size_t i = 1; i < points.size(); i++) {
            const Point& p = points[indices[i]];
            if (!c.contains(p)) {
                c = Circle(p, 0);
                for (size_t j = 0; j < i; j++) {
                    const Point& q = points[indices[j]];
                    if (!c.contains(q)) {
                        c = Circle(p, q);
                        for (size_t k = 0; k < j; k++) {
                            const Point& r = points[indices[k]];
                            if (!c.contains(r)) {
                                c = Circle(p, q, r);
                            }
                        }
                    }
                }
            }
        }

        return c;
    }
};

class DelaunayTriangulation {
private:
    vector<Point> points;
    vector<Triangle> triangles;

    bool isCCW(const Point& a, const Point& b, const Point& c) const {
        return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x) > EPS;
    }

    void addTriangle(const Point& a, const Point& b, const Point& c, int ia, int ib, int ic) {
        if (isCCW(a, b, c)) {
            triangles.push_back(Triangle(a, b, c, ia, ib, ic));
        } else {
            triangles.push_back(Triangle(a, c, b, ia, ic, ib));
        }
    }

    void fixIndices() {
        for (size_t t = 0; t < triangles.size(); t++) {
            Triangle& tri = triangles[t];
            if (tri.ia < 0) {
                for (size_t i = 0; i < points.size(); i++) {
                    if (points[i] == tri.a) { tri.ia = i; break; }
                }
            }
            if (tri.ib < 0) {
                for (size_t i = 0; i < points.size(); i++) {
                    if (points[i] == tri.b) { tri.ib = i; break; }
                }
            }
            if (tri.ic < 0) {
                for (size_t i = 0; i < points.size(); i++) {
                    if (points[i] == tri.c) { tri.ic = i; break; }
                }
            }
        }
    }

public:
    DelaunayTriangulation(const vector<Point>& pts) : points(pts) {}

    void triangulate() {
        if (points.size() < 3) return;

        vector<int> indices(points.size());
        for (size_t i = 0; i < points.size(); i++) indices[i] = i;
        sort(indices.begin(), indices.end(),
             [&](int i, int j) { return points[i] < points[j]; });

        vector<Point> sortedPts;
        vector<int> sortedIndices;
        for (size_t i = 0; i < indices.size(); i++) {
            sortedPts.push_back(points[indices[i]]);
            sortedIndices.push_back(indices[i]);
        }

        double maxX = sortedPts[0].x, maxY = sortedPts[0].y;
        double minX = sortedPts[0].x, minY = sortedPts[0].y;
        for (size_t i = 0; i < sortedPts.size(); i++) {
            const Point& p = sortedPts[i];
            maxX = max(maxX, p.x); maxY = max(maxY, p.y);
            minX = min(minX, p.x); minY = min(minY, p.y);
        }
        double dx = maxX - minX + 1;
        double dy = maxY - minY + 1;

        Point superA(minX - dx * 2, minY - dy * 2);
        Point superB(maxX + dx * 2, minY - dy * 2);
        Point superC((minX + maxX) / 2, maxY + dy * 2);

        triangles.push_back(Triangle(superA, superB, superC, -1, -1, -1));

        for (size_t i = 0; i < sortedPts.size(); i++) {
            Point p = sortedPts[i];
            int originalIdx = sortedIndices[i];

            vector<Triangle> badTriangles;
            for (size_t j = 0; j < triangles.size(); j++) {
                if (triangles[j].containsInCircumcircle(p)) {
                    badTriangles.push_back(triangles[j]);
                }
            }

            vector<pair<Point, Point> > boundary;
            for (size_t j = 0; j < badTriangles.size(); j++) {
                vector<pair<Point, Point> > edges;
                edges.push_back(make_pair(badTriangles[j].a, badTriangles[j].b));
                edges.push_back(make_pair(badTriangles[j].b, badTriangles[j].c));
                edges.push_back(make_pair(badTriangles[j].c, badTriangles[j].a));

                for (size_t e = 0; e < edges.size(); e++) {
                    pair<Point, Point>& edge = edges[e];
                    bool shared = false;
                    for (size_t k = 0; k < badTriangles.size(); k++) {
                        if (j == k) continue;
                        vector<pair<Point, Point> > otherEdges;
                        otherEdges.push_back(make_pair(badTriangles[k].a, badTriangles[k].b));
                        otherEdges.push_back(make_pair(badTriangles[k].b, badTriangles[k].c));
                        otherEdges.push_back(make_pair(badTriangles[k].c, badTriangles[k].a));

                        for (size_t oe = 0; oe < otherEdges.size(); oe++) {
                            pair<Point, Point>& other = otherEdges[oe];
                            if ((edge.first == other.first && edge.second == other.second) ||
                                (edge.first == other.second && edge.second == other.first)) {
                                shared = true;
                                break;
                            }
                        }
                        if (shared) break;
                    }
                    if (!shared) boundary.push_back(edge);
                }
            }

            vector<Triangle> newTriangles;
            for (size_t j = 0; j < triangles.size(); j++) {
                bool isBad = false;
                for (size_t k = 0; k < badTriangles.size(); k++) {
                    if (triangles[j].a == badTriangles[k].a &&
                        triangles[j].b == badTriangles[k].b &&
                        triangles[j].c == badTriangles[k].c) {
                        isBad = true;
                        break;
                    }
                }
                if (!isBad) newTriangles.push_back(triangles[j]);
            }
            triangles = newTriangles;

            for (size_t j = 0; j < boundary.size(); j++) {
                pair<Point, Point>& edge = boundary[j];
                int idx1 = -1, idx2 = -1;
                for (size_t k = 0; k < points.size(); k++) {
                    if (points[k] == edge.first) idx1 = k;
                    if (points[k] == edge.second) idx2 = k;
                }
                addTriangle(edge.first, edge.second, p, idx1, idx2, originalIdx);
            }
        }

        vector<Triangle> filtered;
        for (size_t i = 0; i < triangles.size(); i++) {
            const Triangle& tri = triangles[i];
            bool hasSuper = false;
            if ((fabs(tri.a.x - superA.x) < EPS && fabs(tri.a.y - superA.y) < EPS) ||
                (fabs(tri.b.x - superA.x) < EPS && fabs(tri.b.y - superA.y) < EPS) ||
                (fabs(tri.c.x - superA.x) < EPS && fabs(tri.c.y - superA.y) < EPS) ||
                (fabs(tri.a.x - superB.x) < EPS && fabs(tri.a.y - superB.y) < EPS) ||
                (fabs(tri.b.x - superB.x) < EPS && fabs(tri.b.y - superB.y) < EPS) ||
                (fabs(tri.c.x - superB.x) < EPS && fabs(tri.c.y - superB.y) < EPS) ||
                (fabs(tri.a.x - superC.x) < EPS && fabs(tri.a.y - superC.y) < EPS) ||
                (fabs(tri.b.x - superC.x) < EPS && fabs(tri.b.y - superC.y) < EPS) ||
                (fabs(tri.c.x - superC.x) < EPS && fabs(tri.c.y - superC.y) < EPS)) {
                hasSuper = true;
            }
            if (!hasSuper) filtered.push_back(tri);
        }
        triangles = filtered;

        fixIndices();
    }

    const vector<Triangle>& getTriangles() const {
        return triangles;
    }

    const vector<Point>& getPoints() const {
        return points;
    }
};

DCEL buildVoronoiFromDelaunay(const vector<Triangle>& triangles) {
    DCEL dcel;

    if (triangles.empty()) return dcel;
    map<pair<int, int>, int> edgeToTriangle;
    map<pair<int, int>, int> edgeCount;

    for (size_t i = 0; i < triangles.size(); i++) {
        const Triangle& tri = triangles[i];
        vector<pair<int, int> > edges;
        edges.push_back(make_pair(tri.ia, tri.ib));
        edges.push_back(make_pair(tri.ib, tri.ic));
        edges.push_back(make_pair(tri.ic, tri.ia));

        for (size_t e = 0; e < edges.size(); e++) {
            const pair<int, int>& edge = edges[e];
            if (edge.first < 0 || edge.second < 0) continue;
            int a = min(edge.first, edge.second);
            int b = max(edge.first, edge.second);
            pair<int, int> key = make_pair(a, b);
            edgeToTriangle[key] = i;
            edgeCount[key]++;
        }
    }

    set<pair<int, int> > processedEdges;

    for (size_t i = 0; i < triangles.size(); i++) {
        const Triangle& tri = triangles[i];
        vector<pair<Point, Point> > pointEdges;
        pointEdges.push_back(make_pair(tri.a, tri.b));
        pointEdges.push_back(make_pair(tri.b, tri.c));
        pointEdges.push_back(make_pair(tri.c, tri.a));

        vector<pair<int, int> > indexEdges;
        indexEdges.push_back(make_pair(tri.ia, tri.ib));
        indexEdges.push_back(make_pair(tri.ib, tri.ic));
        indexEdges.push_back(make_pair(tri.ic, tri.ia));

        for (int e = 0; e < 3; e++) {
            const pair<int, int>& indexEdge = indexEdges[e];

            if (indexEdge.first < 0 || indexEdge.second < 0) continue;

            int a = min(indexEdge.first, indexEdge.second);
            int b = max(indexEdge.first, indexEdge.second);
            pair<int, int> key = make_pair(a, b);

            if (edgeCount[key] == 2) {
                int neighborIdx = -1;
                map<pair<int, int>, int>::iterator it = edgeToTriangle.find(key);
                if (it != edgeToTriangle.end() && it->second != (int)i) {
                    neighborIdx = it->second;
                }

                if (neighborIdx != -1) {
                    int id1 = min((int)i, neighborIdx);
                    int id2 = max((int)i, neighborIdx);
                    if (processedEdges.count(make_pair(id1, id2))) continue;
                    processedEdges.insert(make_pair(id1, id2));

                    Point p1 = triangles[i].circumcircle.center;
                    Point p2 = triangles[neighborIdx].circumcircle.center;

                    if (p1.dist2(p2) < EPS) continue;

                    Vertex* v1 = dcel.addVertex(p1);
                    Vertex* v2 = dcel.addVertex(p2);
                    dcel.addEdge(v1, v2);
                }
            }
        }
    }

    set<pair<int, int> > externalEdges;

    for (size_t i = 0; i < triangles.size(); i++) {
        const Triangle& tri = triangles[i];
        vector<pair<Point, Point> > pointEdges;
        pointEdges.push_back(make_pair(tri.a, tri.b));
        pointEdges.push_back(make_pair(tri.b, tri.c));
        pointEdges.push_back(make_pair(tri.c, tri.a));

        vector<pair<int, int> > indexEdges;
        indexEdges.push_back(make_pair(tri.ia, tri.ib));
        indexEdges.push_back(make_pair(tri.ib, tri.ic));
        indexEdges.push_back(make_pair(tri.ic, tri.ia));

        for (int e = 0; e < 3; e++) {
            const pair<Point, Point>& pointEdge = pointEdges[e];
            const pair<int, int>& indexEdge = indexEdges[e];

            if (indexEdge.first < 0 || indexEdge.second < 0) continue;

            int a = min(indexEdge.first, indexEdge.second);
            int b = max(indexEdge.first, indexEdge.second);
            pair<int, int> key = make_pair(a, b);

            if (edgeCount[key] == 1) {
                if (externalEdges.count(key)) continue;
                externalEdges.insert(key);

                Point center = triangles[i].circumcircle.center;

                Point midPoint = Point((pointEdge.first.x + pointEdge.second.x) / 2,
                                      (pointEdge.first.y + pointEdge.second.y) / 2);

                Point dir = midPoint - center;
                double len = sqrt(dir.x*dir.x + dir.y*dir.y);
                if (len > EPS) {
                    dir.x /= len;
                    dir.y /= len;
                } else {
                    Point edgeDir = pointEdge.second - pointEdge.first;
                    dir = Point(-edgeDir.y, edgeDir.x);
                    len = sqrt(dir.x*dir.x + dir.y*dir.y);
                    if (len > EPS) {
                        dir.x /= len;
                        dir.y /= len;
                    } else {
                        dir = Point(1, 0);
                    }
                }

                double INF = 1000.0;
                Point farPoint = center + Point(dir.x * INF, dir.y * INF);

                Vertex* v1 = dcel.addVertex(center);
                Vertex* v2 = dcel.addVertex(farPoint);
                dcel.addEdge(v1, v2);
            }
        }
    }

    return dcel;
}

void printDelaunayTriangles(const vector<Triangle>& triangles) {
    cout << "Треугольники Делоне:" << endl;
    for (size_t i = 0; i < triangles.size(); i++) {
        const Triangle& t = triangles[i];
        cout << "  T" << i << ": (" << t.a.x << ", " << t.a.y << ") ";
        cout << "(" << t.b.x << ", " << t.b.y << ") ";
        cout << "(" << t.c.x << ", " << t.c.y << ")";
        cout << " центр=(" << t.circumcircle.center.x << ", " << t.circumcircle.center.y << ")";
        cout << " R=" << t.circumcircle.radius << endl;
    }
}

void printVoronoiEdges(const DCEL& dcel) {
    cout << "Ребра Вороного:" << endl;

    set<pair<pair<double, double>, pair<double, double> > > uniqueEdges;

    int internalEdges = 0;
    int externalEdges = 0;
    double INF = 500.0;

    for (size_t i = 0; i < dcel.halfEdges.size(); i++) {
        HalfEdge* he = dcel.halfEdges[i];
        if (he->origin && he->twin && he->twin->origin) {
            Point p1 = he->origin->point;
            Point p2 = he->twin->origin->point;

            if (p1.dist2(p2) < EPS) continue;

            pair<double, double> pp1, pp2;
            if (p1.x < p2.x - EPS || (fabs(p1.x - p2.x) < EPS && p1.y < p2.y - EPS)) {
                pp1 = make_pair(p1.x, p1.y);
                pp2 = make_pair(p2.x, p2.y);
            } else {
                pp1 = make_pair(p2.x, p2.y);
                pp2 = make_pair(p1.x, p1.y);
            }

            if (uniqueEdges.count(make_pair(pp1, pp2))) continue;
            uniqueEdges.insert(make_pair(pp1, pp2));

            bool isExternal = (fabs(p1.x) > INF || fabs(p1.y) > INF ||
                              fabs(p2.x) > INF || fabs(p2.y) > INF);

            cout << "  E" << uniqueEdges.size() - 1 << ": (" << p1.x << ", " << p1.y << ") -> (";
            cout << p2.x << ", " << p2.y << ")";
            if (isExternal) {
                cout << " [луч в бесконечность]";
                externalEdges++;
            } else {
                cout << " [внутреннее]";
                internalEdges++;
            }
            cout << endl;
        }
    }

    if (uniqueEdges.empty()) {
        cout << "  Нет ребер" << endl;
    }
    cout << endl << "Статистика:" << endl;
    cout << "  Всего уникальных ребер: " << uniqueEdges.size() << endl;
    cout << "  Внутренних ребер: " << internalEdges << endl;
    cout << "  Лучей в бесконечность: " << externalEdges << endl;
    cout << "  Вершин: " << dcel.vertices.size() << endl;
    cout << "  Полуребер: " << dcel.halfEdges.size() << endl;
}

void visualize(const vector<Point>& points, const vector<Triangle>& triangles, const DCEL& dcel) {
    if (triangles.empty() || points.empty()) {
        cerr << "Нет данных для визуализации" << endl;
        return;
    }

    double minX = points[0].x, maxX = points[0].x;
    double minY = points[0].y, maxY = points[0].y;
    for (size_t i = 0; i < points.size(); i++) {
        const Point& p = points[i];
        minX = min(minX, p.x);
        maxX = max(maxX, p.x);
        minY = min(minY, p.y);
        maxY = max(maxY, p.y);
    }

    double rangeX = maxX - minX;
    double rangeY = maxY - minY;
    double margin = max(rangeX, rangeY) * 0.3 + 2.0;
    if (rangeX < EPS && rangeY < EPS) margin = 5.0;

    minX -= margin;
    maxX += margin;
    minY -= margin;
    maxY += margin;

    const int WINDOW_WIDTH = 1000;
    const int WINDOW_HEIGHT = 700;

    auto toScreenX = [&](double x) {
        return (x - minX) / (maxX - minX) * (WINDOW_WIDTH - 40) + 20;
    };
    auto toScreenY = [&](double y) {
        return (maxY - y) / (maxY - minY) * (WINDOW_HEIGHT - 40) + 20;
    };

    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Delaunay Triangulation & Voronoi Diagram");

    sf::Font font;
    bool fontLoaded = false;
    string fontPaths[] = {
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
        "/usr/share/fonts/truetype/arial/arial.ttf",
        "/System/Library/Fonts/Helvetica.ttf",
        "/Windows/Fonts/arial.ttf"
    };

    for (int i = 0; i < 5; i++) {
        if (font.loadFromFile(fontPaths[i])) {
            fontLoaded = true;
            break;
        }
    }

    sf::Text infoText;
    if (fontLoaded) {
        infoText.setFont(font);
        infoText.setCharacterSize(16);
        infoText.setFillColor(sf::Color::Black);
        infoText.setPosition(10, 10);
    }

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
            if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Escape)
                window.close();
        }

        window.clear(sf::Color::White);

        for (size_t i = 0; i < triangles.size(); i++) {
            const Triangle& t = triangles[i];
            Point center = t.circumcircle.center;
            double radius = t.circumcircle.radius;

            if (radius < 1e6 && radius > EPS) {
                double screenRadius = radius * (WINDOW_WIDTH - 40) / (maxX - minX);
                if (screenRadius > 0 && screenRadius < 10000) {
                    sf::CircleShape circle;
                    circle.setRadius(screenRadius);
                    circle.setPosition(toScreenX(center.x) - screenRadius, toScreenY(center.y) - screenRadius);
                    circle.setFillColor(sf::Color::Transparent);
                    circle.setOutlineColor(sf::Color(255, 200, 50, 100));
                    circle.setOutlineThickness(1.5f);
                    window.draw(circle);

                    sf::CircleShape centerDot(3);
                    centerDot.setPosition(toScreenX(center.x) - 3, toScreenY(center.y) - 3);
                    centerDot.setFillColor(sf::Color(255, 200, 50));
                    centerDot.setOutlineColor(sf::Color(200, 150, 0));
                    centerDot.setOutlineThickness(1);
                    window.draw(centerDot);
                }
            }
        }

        for (size_t i = 0; i < triangles.size(); i++) {
            const Triangle& t = triangles[i];
            sf::ConvexShape triangle;
            triangle.setPointCount(3);
            triangle.setPoint(0, sf::Vector2f(toScreenX(t.a.x), toScreenY(t.a.y)));
            triangle.setPoint(1, sf::Vector2f(toScreenX(t.b.x), toScreenY(t.b.y)));
            triangle.setPoint(2, sf::Vector2f(toScreenX(t.c.x), toScreenY(t.c.y)));
            triangle.setFillColor(sf::Color(200, 220, 255, 80));
            triangle.setOutlineColor(sf::Color(50, 100, 200));
            triangle.setOutlineThickness(1.0f);
            window.draw(triangle);
        }

        for (size_t i = 0; i < dcel.halfEdges.size(); i++) {
            HalfEdge* he = dcel.halfEdges[i];
            if (he->origin && he->twin && he->twin->origin) {
                Point p1 = he->origin->point;
                Point p2 = he->twin->origin->point;

                double len = sqrt(p1.dist2(p2));
                double maxLen = 1000.0;

                Point p1_clip = p1;
                Point p2_clip = p2;

                if (len > maxLen) {
                    Point dir = p2 - p1;
                    dir.x /= len;
                    dir.y /= len;
                    p2_clip = p1 + dir * maxLen;
                }

                if (p2_clip.x < minX - 10 || p2_clip.x > maxX + 10 ||
                    p2_clip.y < minY - 10 || p2_clip.y > maxY + 10) {
                    continue;
                }

                sf::Vertex line[] = {
                    sf::Vertex(sf::Vector2f(toScreenX(p1.x), toScreenY(p1.y)), sf::Color(255, 0, 0, 150)),
                    sf::Vertex(sf::Vector2f(toScreenX(p2_clip.x), toScreenY(p2_clip.y)), sf::Color(255, 0, 0, 150))
                };
                window.draw(line, 2, sf::Lines);
            }
        }

        for (size_t i = 0; i < points.size(); i++) {
            const Point& p = points[i];
            sf::CircleShape circle(7);
            circle.setPosition(toScreenX(p.x) - 7, toScreenY(p.y) - 7);
            circle.setFillColor(sf::Color::Red);
            circle.setOutlineColor(sf::Color::Black);
            circle.setOutlineThickness(1.5);
            window.draw(circle);

            if (fontLoaded) {
                sf::Text text;
                text.setFont(font);
                text.setString(to_string(i));
                text.setCharacterSize(14);
                text.setFillColor(sf::Color::Black);
                text.setPosition(toScreenX(p.x) + 10, toScreenY(p.y) - 8);
                window.draw(text);
            }
        }

        if (fontLoaded) {
            string info = "Points: " + to_string(points.size()) +
                         "  Triangles: " + to_string(triangles.size()) +
                         "\nYellow circles - Circumcircles" +
                         "\nRed lines - Voronoi edges" +
                         "\nBlue - Delaunay triangles" +
                         "\nESC - exit";
            infoText.setString(info);
            window.draw(infoText);
        }

        window.display();
    }
}

int main(int argc, char* argv[]) {
    try {
        vector<Point> points;

        if (argc > 1) {
            for (int i = 1; i + 1 < argc; i += 2) {
                if (i + 1 >= argc) {
                    throw runtime_error("Непарное количество аргументов");
                }
                double x, y;
                stringstream ssx(argv[i]), ssy(argv[i + 1]);
                if (!(ssx >> x) || !(ssy >> y)) {
                    throw runtime_error("Некорректный формат координат");
                }
                points.push_back(Point(x, y));
            }
        } else {
            cout << "Генерация случайных точек..." << endl;
            random_device rd;
            mt19937 gen(rd());
            uniform_real_distribution<> dis(-10, 10);
            for (int i = 0; i < 20; i++) {
                points.push_back(Point(dis(gen), dis(gen)));
            }
        }

        if (points.size() < 3) {
            cerr << "Ошибка: нужно минимум 3 точки" << endl;
            cerr << "Использовано: " << argv[0] << " x1 y1 x2 y2 x3 y3 ..." << endl;
            cerr << "или без аргументов для генерации случайных точек" << endl;
            return 1;
        }

        cout << "Входные точки: ";
        for (size_t i = 0; i < points.size(); i++) {
            const Point& p = points[i];
            cout << "(" << p.x << ", " << p.y << ") ";
        }
        cout << endl << endl;

        MinimumEnclosingCircle mec(points);
        Circle minCircle = mec.compute();
        cout << "Минимальная описанная окружность:" << endl;
        cout << "  Центр: (" << minCircle.center.x << ", " << minCircle.center.y << ")" << endl;
        cout << "  Радиус: " << minCircle.radius << endl << endl;

        DelaunayTriangulation delaunay(points);
        delaunay.triangulate();
        vector<Triangle> triangles = delaunay.getTriangles();

        if (triangles.empty()) {
            cout << "Нет треугольников (точки коллинеарны или вырождены)" << endl;
            return 0;
        }

        printDelaunayTriangles(triangles);

        DCEL dcel = buildVoronoiFromDelaunay(triangles);
        printVoronoiEdges(dcel);

        cout << "\nЗапуск визуализации..." << endl;
        visualize(points, triangles, dcel);

    } catch (const exception& e) {
        cerr << "Ошибка: " << e.what() << endl;
        return 1;
    }
    return 0;
}
