#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <set>
#include <queue>
#include <map>
#include <sstream>
#include <SFML/Graphics.hpp>

using namespace std;

const double EPS = 1e-9;

struct Point {
    double x, y;
    Point(double x = 0, double y = 0) : x(x), y(y) {}
    bool operator<(const Point& other) const {
        if (fabs(x - other.x) > EPS) return x < other.x;
        if (fabs(y - other.y) > EPS) return y < other.y;
        return false;
    }
    bool operator==(const Point& other) const {
        return fabs(x - other.x) < EPS && fabs(y - other.y) < EPS;
    }
    Point operator-(const Point& other) const {
        return Point(x - other.x, y - other.y);
    }
    Point operator+(const Point& other) const {
        return Point(x + other.x, y + other.y);
    }
    Point operator*(double s) const {
        return Point(x * s, y * s);
    }
};

struct Edge {
    int v1, v2;
    Edge(int v1 = -1, int v2 = -1) : v1(v1), v2(v2) {
        if (v1 > v2) swap(this->v1, this->v2);
    }
    bool operator<(const Edge& other) const {
        if (v1 != other.v1) return v1 < other.v1;
        return v2 < other.v2;
    }
    bool operator==(const Edge& other) const {
        return v1 == other.v1 && v2 == other.v2;
    }
};

struct Triangle {
    int v1, v2, v3;
    Triangle(int v1 = -1, int v2 = -1, int v3 = -1)
        : v1(v1), v2(v2), v3(v3) {}
    bool hasVertex(int v) const {
        return v1 == v || v2 == v || v3 == v;
    }
    bool hasEdge(int a, int b) const {
        return (v1 == a && v2 == b) || (v1 == b && v2 == a) ||
               (v1 == a && v3 == b) || (v1 == b && v3 == a) ||
               (v2 == a && v3 == b) || (v2 == b && v3 == a);
    }
    vector<Edge> getEdges() const {
        return {Edge(v1, v2), Edge(v2, v3), Edge(v3, v1)};
    }
    int getThirdVertex(int a, int b) const {
        if (hasEdge(a, b)) {
            if (v1 != a && v1 != b) return v1;
            if (v2 != a && v2 != b) return v2;
            if (v3 != a && v3 != b) return v3;
        }
        return -1;
    }
};

class DelaunayTriangulation {
private:
    vector<Point> points;
    vector<Triangle> triangles;
    map<int, vector<int>> adjacency;
    
    double orient(const Point& a, const Point& b, const Point& c) const {
        return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    }
    
    double dist2(const Point& a, const Point& b) const {
        return (a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y);
    }
    
    bool inCircle(const Point& a, const Point& b, const Point& c, const Point& p) const {
        double ax = a.x - p.x, ay = a.y - p.y;
        double bx = b.x - p.x, by = b.y - p.y;
        double cx = c.x - p.x, cy = c.y - p.y;
        double det = (ax*ax + ay*ay) * (bx*cy - by*cx) -
                     (bx*bx + by*by) * (ax*cy - ay*cx) +
                     (cx*cx + cy*cy) * (ax*by - ay*bx);
        double o = orient(a, b, c);
        if (o < 0) det = -det;
        return det > EPS;
    }

    bool pointInTriangle(const Point& p, const Triangle& t) const {
        Point a = points[t.v1], b = points[t.v2], c = points[t.v3];
        double o1 = orient(a, b, p);
        double o2 = orient(b, c, p);
        double o3 = orient(c, a, p);
        bool hasNeg = (o1 < -EPS) || (o2 < -EPS) || (o3 < -EPS);
        bool hasPos = (o1 > EPS) || (o2 > EPS) || (o3 > EPS);
        return !(hasNeg && hasPos);
    }

    int findTriangleContainingPoint(const Point& p) const {
        for (size_t i = 0; i < triangles.size(); ++i) {
            if (pointInTriangle(p, triangles[i])) {
                return i;
            }
        }
        double minDist = 1e100;
        int bestIdx = 0;
        for (size_t i = 0; i < triangles.size(); ++i) {
            const Triangle& t = triangles[i];
            Point center(
                (points[t.v1].x + points[t.v2].x + points[t.v3].x) / 3.0,
                (points[t.v1].y + points[t.v2].y + points[t.v3].y) / 3.0
            );
            double d = dist2(center, p);
            if (d < minDist) {
                minDist = d;
                bestIdx = i;
            }
        }
        return bestIdx;
    }

    vector<Triangle> splitTriangle(const Triangle& t, int newVertex) {
        return {
            Triangle(t.v1, t.v2, newVertex),
            Triangle(t.v2, t.v3, newVertex),
            Triangle(t.v3, t.v1, newVertex)
        };
    }

    bool isConvexQuadrilateral(int v1, int v2, int v3, int v4) const {
        double o1 = orient(points[v1], points[v2], points[v3]);
        double o2 = orient(points[v1], points[v2], points[v4]);
        return o1 * o2 < -EPS;
    }

    void buildAdjacency() {
        adjacency.clear();
        for (int i = 0; i < (int)points.size(); ++i) {
            adjacency[i] = vector<int>();
        }
        for (const Triangle& t : triangles) {
            adjacency[t.v1].push_back(t.v2);
            adjacency[t.v1].push_back(t.v3);
            adjacency[t.v2].push_back(t.v1);
            adjacency[t.v2].push_back(t.v3);
            adjacency[t.v3].push_back(t.v1);
            adjacency[t.v3].push_back(t.v2);
        }
    }

public:
    DelaunayTriangulation() {}

    const vector<Triangle>& getTriangles() const {
        return triangles;
    }
    
    const vector<Point>& getPoints() const {
        return points;
    }

    void addPoints(const vector<Point>& newPoints) {
        points = newPoints;
        buildTriangulation();
    }

    void buildTriangulation() {
        if (points.size() < 3) {
            triangles.clear();
            return;
        }
      
        int originalSize = points.size();
        double minX = points[0].x, maxX = points[0].x;
        double minY = points[0].y, maxY = points[0].y;
        for (const Point& p : points) {
            minX = min(minX, p.x);
            maxX = max(maxX, p.x);
            minY = min(minY, p.y);
            maxY = max(maxY, p.y);
        }
      
        double dx = maxX - minX;
        double dy = maxY - minY;
        double d = max(dx, dy) * 10.0 + 100.0;
        double cx = (minX + maxX) / 2.0;
        double cy = (minY + maxY) / 2.0;
        Point p1(cx - d * 2, cy - d);
        Point p2(cx + d * 2, cy - d);
        Point p3(cx, cy + d * 3);
        points.push_back(p1);
        points.push_back(p2);
        points.push_back(p3);
        int superV1 = originalSize;
        int superV2 = originalSize + 1;
        int superV3 = originalSize + 2;
        triangles.clear();
        triangles.push_back(Triangle(superV1, superV2, superV3));

        for (int i = 0; i < originalSize; ++i) {
            insertPoint(i);
        }

        vector<Triangle> filtered;
        for (const Triangle& t : triangles) {
            bool hasSuper = t.hasVertex(superV1) || t.hasVertex(superV2) || t.hasVertex(superV3);
            if (!hasSuper) {
                filtered.push_back(t);
            }
        }
        triangles = filtered;
        points.resize(originalSize);

        if (triangles.empty() && points.size() >= 3) {
            if (points.size() == 3) {
                triangles.push_back(Triangle(0, 1, 2));
            } else {
                vector<int> indices(points.size());
                for (size_t i = 0; i < points.size(); ++i) indices[i] = i;
                double cx_sum = 0, cy_sum = 0;
                for (const Point& p : points) {
                    cx_sum += p.x;
                    cy_sum += p.y;
                }
                cx_sum /= points.size();
                cy_sum /= points.size();
                sort(indices.begin(), indices.end(), [&](int a, int b) {
                    double angle_a = atan2(points[a].y - cy_sum, points[a].x - cx_sum);
                    double angle_b = atan2(points[b].y - cy_sum, points[b].x - cx_sum);
                    return angle_a < angle_b;
                });
                for (int i = 1; i < (int)indices.size() - 1; ++i) {
                    triangles.push_back(Triangle(indices[0], indices[i], indices[i + 1]));
                }
            }
        }
    }

    void insertPoint(int pointIdx) {
        Point p = points[pointIdx];
        int triIdx = findTriangleContainingPoint(p);

        if (triIdx == -1 || triIdx >= (int)triangles.size()) {
            triIdx = 0;
        }

        Triangle t = triangles[triIdx];
        triangles.erase(triangles.begin() + triIdx);
        vector<Triangle> newTriangles = splitTriangle(t, pointIdx);
        triangles.insert(triangles.end(), newTriangles.begin(), newTriangles.end());
        set<Edge> edgesToCheck;
        for (const Triangle& nt : newTriangles) {
            vector<Edge> edges = nt.getEdges();
            for (const Edge& e : edges) {
                edgesToCheck.insert(e);
            }
        }

        queue<Edge> edgeQueue;
        for (const Edge& e : edgesToCheck) {
            edgeQueue.push(e);
        }

        set<Edge> processedEdges;
        int maxIterations = 10000;
        int iterations = 0;

        while (!edgeQueue.empty() && iterations < maxIterations) {
            Edge e = edgeQueue.front();
            edgeQueue.pop();
            iterations++;

            if (processedEdges.count(e)) continue;
            processedEdges.insert(e);

            int tri1 = -1, tri2 = -1;
            int v3 = -1, v4 = -1;

            for (size_t i = 0; i < triangles.size(); ++i) {
                if (triangles[i].hasEdge(e.v1, e.v2)) {
                    if (tri1 == -1) {
                        tri1 = i;
                        v3 = triangles[i].getThirdVertex(e.v1, e.v2);
                    } else {
                        tri2 = i;
                        v4 = triangles[i].getThirdVertex(e.v1, e.v2);
                        break;
                    }
                }
            }

            if (tri1 != -1 && tri2 != -1 && v3 != -1 && v4 != -1) {
                if (inCircle(points[e.v1], points[e.v2], points[v3], points[v4]) ||
                    inCircle(points[e.v1], points[e.v2], points[v4], points[v3])) {

                    if (isConvexQuadrilateral(e.v1, e.v2, v3, v4)) {
                        triangles[tri1] = Triangle(e.v1, v3, v4);
                        triangles[tri2] = Triangle(e.v2, v3, v4);

                        Edge e1(e.v1, v3);
                        Edge e2(e.v1, v4);
                        Edge e3(e.v2, v3);
                        Edge e4(e.v2, v4);

                        if (!processedEdges.count(e1)) edgeQueue.push(e1);
                        if (!processedEdges.count(e2)) edgeQueue.push(e2);
                        if (!processedEdges.count(e3)) edgeQueue.push(e3);
                        if (!processedEdges.count(e4)) edgeQueue.push(e4);
                    }
                }
            }
        }
    }

    void validatePoints(const vector<Point>& pts) const {
        if (pts.size() < 3) {
            throw runtime_error("Для триангуляции нужно минимум 3 точки");
        }
        for (size_t i = 0; i < pts.size(); ++i) {
            for (size_t j = i + 1; j < pts.size(); ++j) {
                if (pts[i] == pts[j]) {
                    throw runtime_error("Обнаружены дублирующиеся точки");
                }
            }
        }
        bool allCollinear = true;
        for (size_t i = 2; i < pts.size(); ++i) {
            if (fabs(orient(pts[0], pts[1], pts[i])) > EPS) {
                allCollinear = false;
                break;
            }
        }
        if (allCollinear && pts.size() > 2) {
            throw runtime_error("Все точки лежат на одной прямой");
        }
    }

    void printTriangulation() const {
        if (triangles.empty()) {
            cout << "Триангуляция не выполнена" << endl;
            return;
        }
        cout << "Треугольники Делоне:" << endl;
        for (size_t i = 0; i < triangles.size(); ++i) {
            const Triangle& t = triangles[i];
            cout << "  T" << i << ": (" << points[t.v1].x << ", " << points[t.v1].y << ") ";
            cout << "(" << points[t.v2].x << ", " << points[t.v2].y << ") ";
            cout << "(" << points[t.v3].x << ", " << points[t.v3].y << ")";
            cout << " индексы: (" << t.v1 << ", " << t.v2 << ", " << t.v3 << ")" << endl;
        }
        cout << "Всего треугольников: " << triangles.size() << endl;
    }

    void visualize() const {
        if (triangles.empty()) {
            cerr << "Триангуляция пуста, визуализация невозможна" << endl;
            return;
        }

        double minX = points[0].x, maxX = points[0].x;
        double minY = points[0].y, maxY = points[0].y;
        for (const Point& p : points) {
            minX = min(minX, p.x);
            maxX = max(maxX, p.x);
            minY = min(minY, p.y);
            maxY = max(maxY, p.y);
        }

        double margin = max(maxX - minX, maxY - minY) * 0.2 + 1.0;
        minX -= margin;
        maxX += margin;
        minY -= margin;
        maxY += margin;
        const int WINDOW_WIDTH = 800;
        const int WINDOW_HEIGHT = 600;

        auto toScreenX = [&](double x) {
            return (x - minX) / (maxX - minX) * (WINDOW_WIDTH - 40) + 20;
        };
        auto toScreenY = [&](double y) {
            return (maxY - y) / (maxY - minY) * (WINDOW_HEIGHT - 40) + 20;
        };

        sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Триангуляция Делоне");
        
        sf::Font font;
        if (!font.loadFromFile("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf")) {
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

            for (const Triangle& t : triangles) {
                Point p1 = points[t.v1];
                Point p2 = points[t.v2];
                Point p3 = points[t.v3];

                sf::ConvexShape triangle;
                triangle.setPointCount(3);
                triangle.setPoint(0, sf::Vector2f(toScreenX(p1.x), toScreenY(p1.y)));
                triangle.setPoint(1, sf::Vector2f(toScreenX(p2.x), toScreenY(p2.y)));
                triangle.setPoint(2, sf::Vector2f(toScreenX(p3.x), toScreenY(p3.y)));
                triangle.setFillColor(sf::Color(200, 220, 255, 100));
                triangle.setOutlineColor(sf::Color(50, 100, 200));
                triangle.setOutlineThickness(1.5f);
                window.draw(triangle);
            }

            for (size_t i = 0; i < points.size(); ++i) {
                const Point& p = points[i];
                sf::CircleShape circle(6);
                circle.setPosition(toScreenX(p.x) - 6, toScreenY(p.y) - 6);
                circle.setFillColor(sf::Color::Red);
                circle.setOutlineColor(sf::Color::Black);
                circle.setOutlineThickness(1);
                window.draw(circle);

                sf::Text text;
                text.setFont(font);
                text.setString(to_string(i));
                text.setCharacterSize(14);
                text.setFillColor(sf::Color::Black);
                text.setPosition(toScreenX(p.x) + 8, toScreenY(p.y) - 8);
                window.draw(text);
            }

            sf::Text info;
            info.setFont(font);
            info.setString("Points: " + to_string(points.size()) + 
                          "  Triangles: " + to_string(triangles.size()) +
                          "\nESC - quit");
            info.setCharacterSize(16);
            info.setFillColor(sf::Color::Black);
            info.setPosition(10, 10);
            window.draw(info);

            window.display();
        }
    }
};

vector<Point> readPointsFromArgs(int argc, char* argv[]) {
    vector<Point> points;
    for (int i = 1; i < argc; i += 2) {
        if (i + 1 >= argc) {
            throw runtime_error("Непарное количество аргументов");
        }
        double x, y;
        stringstream ssx(argv[i]), ssy(argv[i + 1]);
        if (!(ssx >> x) || !(ssy >> y)) {
            throw runtime_error("Неправильный формат координат: " + string(argv[i]) + " " + string(argv[i+1]));
        }
        points.push_back(Point(x, y));
    }
    return points;
}

int main(int argc, char* argv[]) {
    try {
        vector<Point> points = readPointsFromArgs(argc, argv);

        cout << "Входные точки: ";
        for (const auto& p : points) {
            cout << "(" << p.x << ", " << p.y << ") ";
        }
        cout << endl << endl;

        DelaunayTriangulation dt;
        dt.validatePoints(points);
        dt.addPoints(points);
        dt.printTriangulation();
        dt.visualize();
        
    } catch (const exception& e) {
        cerr << "Ошибка: " << e.what() << endl;
        return 1;
    }
    return 0;
}
