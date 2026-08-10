#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <set>
#include <sstream>
#include <limits>

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
        for (auto v : vertices) delete v;
        for (auto he : halfEdges) delete he;
        for (auto f : faces) delete f;
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
        HalfEdge* he1 = addHalfEdge();
        HalfEdge* he2 = addHalfEdge();
        he1->origin = v1;
        he2->origin = v2;
        he1->twin = he2;
        he2->twin = he1;
        if (!v1->incident) v1->incident = he1;
        if (!v2->incident) v2->incident = he2;
    }
};

struct Triangle {
    Point a, b, c;
    int ia, ib, ic;
    Point circumcenter;

    Triangle(Point a, Point b, Point c, int ia, int ib, int ic)
        : a(a), b(b), c(c), ia(ia), ib(ib), ic(ic) {
        updateCircumcenter();
    }

    void updateCircumcenter() {
        double d = 2 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
        if (fabs(d) < EPS) {
            circumcenter = Point(0, 0);
            return;
        }
        double ux = (a.x*a.x + a.y*a.y) * (b.y - c.y) +
                    (b.x*b.x + b.y*b.y) * (c.y - a.y) +
                    (c.x*c.x + c.y*c.y) * (a.y - b.y);
        double uy = (a.x*a.x + a.y*a.y) * (c.x - b.x) +
                    (b.x*b.x + b.y*b.y) * (a.x - c.x) +
                    (c.x*c.x + c.y*c.y) * (b.x - a.x);
        circumcenter = Point(ux / d, uy / d);
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
};

class DelaunayTriangulation {
private:
    vector<Point> points;
    vector<Triangle> triangles;

    int findPointIndex(const Point& p) {
        for (int i = 0; i < points.size(); i++) {
            if (points[i] == p) return i;
        }
        return -1;
    }

    bool isCCW(const Point& a, const Point& b, const Point& c) {
        return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x) > EPS;
    }

    void addTriangle(const Point& a, const Point& b, const Point& c, int ia, int ib, int ic) {
        if (isCCW(a, b, c)) {
            triangles.push_back(Triangle(a, b, c, ia, ib, ic));
        } else {
            triangles.push_back(Triangle(a, c, b, ia, ic, ib));
        }
    }

public:
    DelaunayTriangulation(const vector<Point>& pts) : points(pts) {}

    void triangulate() {
        if (points.size() < 3) return;

        vector<int> indices(points.size());
        for (int i = 0; i < points.size(); i++) indices[i] = i;
        sort(indices.begin(), indices.end(),
             [&](int i, int j) { return points[i] < points[j]; });

        vector<Point> sortedPts;
        vector<int> sortedIndices;
        for (int i : indices) {
            sortedPts.push_back(points[i]);
            sortedIndices.push_back(i);
        }

        double maxX = sortedPts[0].x, maxY = sortedPts[0].y;
        double minX = sortedPts[0].x, minY = sortedPts[0].y;
        for (const auto& p : sortedPts) {
            maxX = max(maxX, p.x); maxY = max(maxY, p.y);
            minX = min(minX, p.x); minY = min(minY, p.y);
        }
        double dx = maxX - minX + 1;
        double dy = maxY - minY + 1;

        Point superA(minX - dx * 2, minY - dy * 2);
        Point superB(maxX + dx * 2, minY - dy * 2);
        Point superC((minX + maxX) / 2, maxY + dy * 2);

        triangles.push_back(Triangle(superA, superB, superC, -1, -1, -1));

        for (int i = 0; i < sortedPts.size(); i++) {
            Point p = sortedPts[i];
            int originalIdx = sortedIndices[i];

            vector<Triangle> badTriangles;
            for (const auto& tri : triangles) {
                Point center = tri.circumcenter;
                if (fabs(center.x) < EPS && fabs(center.y) < EPS) {
                    badTriangles.push_back(tri);
                    continue;
                }
                double r2 = tri.a.dist2(center);
                double d2 = p.dist2(center);
                if (d2 < r2 - EPS) {
                    badTriangles.push_back(tri);
                }
            }

            vector<pair<Point, Point>> boundary;
            for (int j = 0; j < badTriangles.size(); j++) {
                vector<pair<Point, Point>> edges = {
                    {badTriangles[j].a, badTriangles[j].b},
                    {badTriangles[j].b, badTriangles[j].c},
                    {badTriangles[j].c, badTriangles[j].a}
                };
                for (auto& e : edges) {
                    bool shared = false;
                    for (int k = 0; k < badTriangles.size(); k++) {
                        if (j == k) continue;
                        vector<pair<Point, Point>> otherEdges = {
                            {badTriangles[k].a, badTriangles[k].b},
                            {badTriangles[k].b, badTriangles[k].c},
                            {badTriangles[k].c, badTriangles[k].a}
                        };
                        for (auto& oe : otherEdges) {
                            if ((e.first == oe.first && e.second == oe.second) ||
                                (e.first == oe.second && e.second == oe.first)) {
                                shared = true;
                                break;
                            }
                        }
                        if (shared) break;
                    }
                    if (!shared) boundary.push_back(e);
                }
            }

            vector<Triangle> newTriangles;
            for (const auto& tri : triangles) {
                bool isBad = false;
                for (const auto& bad : badTriangles) {
                    if (tri.a == bad.a && tri.b == bad.b && tri.c == bad.c) {
                        isBad = true;
                        break;
                    }
                }
                if (!isBad) newTriangles.push_back(tri);
            }
            triangles = newTriangles;

            for (auto& edge : boundary) {
                int idx1 = findPointIndex(edge.first);
                int idx2 = findPointIndex(edge.second);
                if (idx1 == -1) idx1 = -1;
                if (idx2 == -1) idx2 = -1;
                addTriangle(edge.first, edge.second, p, idx1, idx2, originalIdx);
            }
        }

        vector<Triangle> filtered;
        for (const auto& tri : triangles) {
            bool hasSuper = false;
            for (const auto& v : {tri.a, tri.b, tri.c}) {
                if ((fabs(v.x - superA.x) < EPS && fabs(v.y - superA.y) < EPS) ||
                    (fabs(v.x - superB.x) < EPS && fabs(v.y - superB.y) < EPS) ||
                    (fabs(v.x - superC.x) < EPS && fabs(v.y - superC.y) < EPS)) {
                    hasSuper = true;
                    break;
                }
            }
            if (!hasSuper) filtered.push_back(tri);
        }
        triangles = filtered;
    }

    const vector<Triangle>& getTriangles() const {
        return triangles;
    }
};

DCEL buildVoronoiFromDelaunay(const vector<Triangle>& triangles) {
    DCEL dcel;

    if (triangles.empty()) return dcel;

    set<pair<int, int>> processedEdges;
    set<pair<int, int>> externalEdges;

    for (int i = 0; i < triangles.size(); i++) {
        const auto& tri = triangles[i];
        vector<pair<Point, Point>> pointEdges = {
            {tri.a, tri.b}, {tri.b, tri.c}, {tri.c, tri.a}
        };
        vector<pair<int, int>> indexEdges = {
            {tri.ia, tri.ib}, {tri.ib, tri.ic}, {tri.ic, tri.ia}
        };

        for (int e = 0; e < 3; e++) {
            const auto& pointEdge = pointEdges[e];
            const auto& indexEdge = indexEdges[e];

            if (indexEdge.first < 0 || indexEdge.second < 0) continue;

            int neighborIdx = -1;
            for (int j = 0; j < triangles.size(); j++) {
                if (i == j) continue;
                const auto& other = triangles[j];
                if (other.hasEdgePoints(pointEdge.first, pointEdge.second)) {
                    neighborIdx = j;
                    break;
                }
            }

            if (neighborIdx != -1) {
                int id1 = min(i, neighborIdx);
                int id2 = max(i, neighborIdx);
                if (processedEdges.count({id1, id2})) continue;
                processedEdges.insert({id1, id2});

                Point p1 = triangles[i].circumcenter;
                Point p2 = triangles[neighborIdx].circumcenter;

                if (p1.dist2(p2) < EPS) {
                    continue;
                }

                Vertex* v1 = dcel.addVertex(p1);
                Vertex* v2 = dcel.addVertex(p2);
                dcel.addEdge(v1, v2);
            }
        }
    }

    for (int i = 0; i < triangles.size(); i++) {
        const auto& tri = triangles[i];
        vector<pair<Point, Point>> pointEdges = {
            {tri.a, tri.b}, {tri.b, tri.c}, {tri.c, tri.a}
        };
        vector<pair<int, int>> indexEdges = {
            {tri.ia, tri.ib}, {tri.ib, tri.ic}, {tri.ic, tri.ia}
        };

        for (int e = 0; e < 3; e++) {
            const auto& pointEdge = pointEdges[e];
            const auto& indexEdge = indexEdges[e];

            if (indexEdge.first < 0 || indexEdge.second < 0) continue;

            bool hasNeighbor = false;
            for (int j = 0; j < triangles.size(); j++) {
                if (i == j) continue;
                const auto& other = triangles[j];
                if (other.hasEdgePoints(pointEdge.first, pointEdge.second)) {
                    hasNeighbor = true;
                    break;
                }
            }

            if (!hasNeighbor) {
                int id1 = min(indexEdge.first, indexEdge.second);
                int id2 = max(indexEdge.first, indexEdge.second);
                if (externalEdges.count({id1, id2})) continue;
                externalEdges.insert({id1, id2});

                Point center = triangles[i].circumcenter;

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
    for (int i = 0; i < triangles.size(); i++) {
        const auto& t = triangles[i];
        cout << "  T" << i << ": (" << t.a.x << ", " << t.a.y << ") ";
        cout << "(" << t.b.x << ", " << t.b.y << ") ";
        cout << "(" << t.c.x << ", " << t.c.y << ")";
        cout << " индексы: (" << t.ia << ", " << t.ib << ", " << t.ic << ")" << endl;
    }
}

void printVoronoiEdges(const DCEL& dcel) {
    cout << "Ребра Вороного:" << endl;

    set<pair<pair<double, double>, pair<double, double>>> uniqueEdges;

    int internalEdges = 0;
    int externalEdges = 0;
    double INF = 500.0;

    for (auto he : dcel.halfEdges) {
        if (he->origin && he->twin && he->twin->origin) {
            Point p1 = he->origin->point;
            Point p2 = he->twin->origin->point;

            if (p1.dist2(p2) < EPS) {
                continue;
            }

            pair<double, double> pp1, pp2;
            if (p1.x < p2.x - EPS || (fabs(p1.x - p2.x) < EPS && p1.y < p2.y - EPS)) {
                pp1 = {p1.x, p1.y};
                pp2 = {p2.x, p2.y};
            } else {
                pp1 = {p2.x, p2.y};
                pp2 = {p1.x, p1.y};
            }

            if (uniqueEdges.count({pp1, pp2})) continue;
            uniqueEdges.insert({pp1, pp2});

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
        cout << "  Нет ребер (вырожденный случай)" << endl;
    }
    cout << endl << "Статистика:" << endl;
    cout << "  Всего уникальных ребер: " << uniqueEdges.size() << endl;
    cout << "  Внутренних ребер: " << internalEdges << endl;
    cout << "  Лучей в бесконечность: " << externalEdges << endl;
    cout << "  Вершин: " << dcel.vertices.size() << endl;
    cout << "  Полуребер: " << dcel.halfEdges.size() << endl;
}

int main(int argc, char* argv[]) {
    try {
        vector<Point> points;

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

        if (points.size() < 3) {
            cerr << "Ошибка: нужно минимум 3 точк" << endl;
            cerr << "Использовано: " << argv[0] << " x1 y1 x2 y2 x3 y3 ..." << endl;
            return 1;
        }

        cout << "Входные точки: ";
        for (const auto& p : points) {
            cout << "(" << p.x << ", " << p.y << ") ";
        }
        cout << endl << endl;

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
    } catch (const exception& e) {
        cerr << "Ошибка: " << e.what() << endl;
        return 1;
    }
    return 0;
}
