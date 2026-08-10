#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <set>
#include <stack>
#include <sstream>

using namespace std;

const double EPS = 1e-9;

struct Point {
    double x, y;
    Point(double x = 0, double y = 0) : x(x), y(y) {}
    Point operator-(const Point& other) const {
        return Point(x - other.x, y - other.y);
    }
    bool operator==(const Point& other) const {
        return fabs(x - other.x) < EPS && fabs(y - other.y) < EPS;
    }
    bool operator<(const Point& other) const {
        if (fabs(x - other.x) > EPS) return x < other.x;
        if (fabs(y - other.y) > EPS) return y < other.y;
        return false;
    }
};

struct Vertex {
    Point p;
    Vertex* prev;
    Vertex* next;
    int index;
    bool removed;
    bool isConvex;
    Vertex(const Point& pt, int idx) : p(pt), index(idx), prev(nullptr), next(nullptr), removed(false), isConvex(false) {}
};

double crossProduct(const Point& a, const Point& b, const Point& c) {
    Point ab = b - a;
    Point ac = c - a;
    return ab.x * ac.y - ab.y * ac.x;
}

bool isCounterClockwise(const vector<Point>& polygon) {
    double area = 0;
    int n = polygon.size();
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        area += polygon[i].x * polygon[j].y;
        area -= polygon[j].x * polygon[i].y;
    }
    return area > EPS;
}

bool isConvex(Vertex* v) {
    return crossProduct(v->prev->p, v->p, v->next->p) >= -EPS;
}

bool pointInTriangle(const Point& p, const Point& a, const Point& b, const Point& c) {
    double d1 = crossProduct(a, b, p);
    double d2 = crossProduct(b, c, p);
    double d3 = crossProduct(c, a, p);
    bool hasNeg = (d1 < -EPS) || (d2 < -EPS) || (d3 < -EPS);
    bool hasPos = (d1 > EPS) || (d2 > EPS) || (d3 > EPS);
    return !(hasNeg && hasPos);
}

bool isEar(Vertex* v) {
    if (!isConvex(v)) return false;
    Vertex* current = v->next->next;
    while (current != v->prev) {
        if (!current->removed && pointInTriangle(current->p, v->prev->p, v->p, v->next->p)) {
            return false;
        }
        current = current->next;
    }
    return true;
}

bool segmentsIntersect(const Point& a, const Point& b, const Point& c, const Point& d) {
    double d1 = crossProduct(a, b, c);
    double d2 = crossProduct(a, b, d);
    double d3 = crossProduct(c, d, a);
    double d4 = crossProduct(c, d, b);
    if (fabs(d1) < EPS && fabs(d2) < EPS && fabs(d3) < EPS && fabs(d4) < EPS) {
        double min1x = min(a.x, b.x), max1x = max(a.x, b.x);
        double min1y = min(a.y, b.y), max1y = max(a.y, b.y);
        double min2x = min(c.x, d.x), max2x = max(c.x, d.x);
        double min2y = min(c.y, d.y), max2y = max(c.y, d.y);

        bool xOverlap = (min1x <= max2x + EPS && max1x >= min2x - EPS);
        bool yOverlap = (min1y <= max2y + EPS && max1y >= min2y - EPS);
        return xOverlap && yOverlap;
    }

    bool intersect1 = (d1 > EPS && d2 < -EPS) || (d1 < -EPS && d2 > EPS);
    bool intersect2 = (d3 > EPS && d4 < -EPS) || (d3 < -EPS && d4 > EPS);

    if (intersect1 && intersect2) return true;

    if (fabs(d1) < EPS && (c.x >= min(a.x, b.x) - EPS && c.x <= max(a.x, b.x) + EPS &&
                           c.y >= min(a.y, b.y) - EPS && c.y <= max(a.y, b.y) + EPS)) return true;
    if (fabs(d2) < EPS && (d.x >= min(a.x, b.x) - EPS && d.x <= max(a.x, b.x) + EPS &&
                           d.y >= min(a.y, b.y) - EPS && d.y <= max(a.y, b.y) + EPS)) return true;
    if (fabs(d3) < EPS && (a.x >= min(c.x, d.x) - EPS && a.x <= max(c.x, d.x) + EPS &&
                           a.y >= min(c.y, d.y) - EPS && a.y <= max(c.y, d.y) + EPS)) return true;
    if (fabs(d4) < EPS && (b.x >= min(c.x, d.x) - EPS && b.x <= max(c.x, d.x) + EPS &&
                           b.y >= min(c.y, d.y) - EPS && b.y <= max(c.y, d.y) + EPS)) return true;

    return false;
}

bool hasSelfIntersections(const vector<Point>& polygon) {
    int n = polygon.size();
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            int next_i = (i + 1) % n;
            int next_j = (j + 1) % n;
            if (j == next_i || i == next_j) continue;
            if (segmentsIntersect(polygon[i], polygon[next_i], polygon[j], polygon[next_j])) {
                return true;
            }
        }
    }
    return false;
}

bool hasDuplicateVertices(const vector<Point>& polygon) {
    for (size_t i = 0; i < polygon.size(); ++i) {
        for (size_t j = i + 1; j < polygon.size(); ++j) {
            if (polygon[i] == polygon[j]) {
                return true;
            }
        }
    }
    return false;
}

bool hasCollinearEdges(const vector<Point>& polygon) {
    int n = polygon.size();
    if (n < 3) return false;
    bool allCollinear = true;
    for (int i = 0; i < n; ++i) {
        Point a = polygon[i];
        Point b = polygon[(i + 1) % n];
        Point c = polygon[(i + 2) % n];
        if (fabs(crossProduct(a, b, c)) > EPS) {
            allCollinear = false;
            break;
        }
    }
    return allCollinear;
}

bool validatePolygon(vector<Point>& polygon, string& errorMsg) {
    int n = polygon.size();
    if (n < 3) {
        errorMsg = "Многоугольник должен содержать как минимум 3 вершины";
        return false;
    }
    if (hasDuplicateVertices(polygon)) {
        errorMsg = "Обнаружены повторяющиеся вершины";
        return false;
    }
    if (hasSelfIntersections(polygon)) {
        errorMsg = "Многоугольник имеет самопересечения";
        return false;
    }
    if (hasCollinearEdges(polygon)) {
        errorMsg = "Все вершины коллинеарны (многоугольник вырожден)";
        return false;
    }
    if (!isCounterClockwise(polygon)) {
        reverse(polygon.begin(), polygon.end());
        cout << "Направление обхода автоматически изменено на против часовой стрелки" << endl;
    }

    return true;
}

vector<Point> triangulate(const vector<Point>& polygon) {
    vector<Point> triangles;
    int n = polygon.size();
    if (n < 3) return triangles;

    vector<Vertex*> vertices(n);
    for (int i = 0; i < n; ++i) {
        vertices[i] = new Vertex(polygon[i], i);
    }

    for (int i = 0; i < n; ++i) {
        vertices[i]->prev = vertices[(i - 1 + n) % n];
        vertices[i]->next = vertices[(i + 1) % n];
    }

    Vertex* head = vertices[0];
    int remaining = n;

    while (remaining > 3) {
        bool earFound = false;
        Vertex* current = head;
        int attempts = 0;

        do {
            if (!current->removed && isEar(current)) {
                Vertex* prev = current->prev;
                Vertex* next = current->next;
                triangles.push_back(prev->p);
                triangles.push_back(current->p);
                triangles.push_back(next->p);
                prev->next = next;
                next->prev = prev;
                current->removed = true;
                if (head == current) {
                    head = next;
                }
                remaining--;
                earFound = true;
                break;
            }
            current = current->next;
            attempts++;
        } while (current != head && attempts < n * 2);

        if (!earFound) {
            current = head;
            attempts = 0;
            do {
                if (!current->removed) {
                    Vertex* prev = current->prev;
                    Vertex* next = current->next;
                    bool hasPoints = false;
                    Vertex* test = next->next;
                    while (test != prev) {
                        if (!test->removed && pointInTriangle(test->p, prev->p, current->p, next->p)) {
                            hasPoints = true;
                            break;
                        }
                        test = test->next;
                    }
                    if (!hasPoints) {
                        triangles.push_back(prev->p);
                        triangles.push_back(current->p);
                        triangles.push_back(next->p);
                        prev->next = next;
                        next->prev = prev;
                        current->removed = true;
                        if (head == current) {
                            head = next;
                        }
                        remaining--;
                        earFound = true;
                        break;
                    }
                }
                current = current->next;
                attempts++;
            } while (current != head && attempts < n * 2);
        }

        if (!earFound) {
            cerr << "Предупреждение: ухо не найдено,завершение" << endl;
            break;
        }
    }

    if (remaining == 3) {
        Vertex* v1 = head;
        Vertex* v2 = head->next;
        Vertex* v3 = head->next->next;
        if (!v1->removed && !v2->removed && !v3->removed) {
            triangles.push_back(v1->p);
            triangles.push_back(v2->p);
            triangles.push_back(v3->p);
        }
    }

    Vertex* current = head;
    for (int i = 0; i < remaining; ++i) {
        Vertex* next = current->next;
        delete current;
        current = next;
    }
    return triangles;
}

vector<Point> readPointsFromArgs(int argc, char* argv[]) {
    vector<Point> points;
    for (int i = 1; i < argc; i += 2) {
        if (i + 1 >= argc) {
            throw runtime_error("Непарное количество аргументов");
        }
        double x, y;
        stringstream ssx(argv[i]), ssy(argv[i + 1]);
        if (!(ssx >> x) || !(ssy >> y)) {
            throw runtime_error("Неправильный формат координат");
        }
        points.push_back(Point(x, y));
    }
    return points;
}

void printTriangles(const vector<Point>& triangles) {
    if (triangles.empty()) {
        cout << "Триангуляция не выполнена" << endl;
        return;
    }

    cout << "Треугольники:" << endl;
    for (size_t i = 0; i < triangles.size(); i += 3) {
        cout << "  T" << i/3 << ": (" << triangles[i].x << ", " << triangles[i].y << ") ";
        cout << "(" << triangles[i+1].x << ", " << triangles[i+1].y << ") ";
        cout << "(" << triangles[i+2].x << ", " << triangles[i+2].y << ")" << endl;
    }
    cout << "Всего треугольников: " << triangles.size() / 3 << endl;
}

int main(int argc, char* argv[]) {
    try {
        if (argc < 7) {
            cout << "Использование: " << argv[0] << " x1 y1 x2 y2 x3 y3 ..." << endl;
            cout << "Пример: " << argv[0] << " 0 0 4 0 4 3 0 3" << endl;
            return 1;
        }

        vector<Point> polygon = readPointsFromArgs(argc, argv);
        string errorMsg;

        cout << "Входной многоугольник: ";
        for (const auto& p : polygon) {
            cout << "(" << p.x << ", " << p.y << ") ";
        }
        cout << endl << endl;

        if (!validatePolygon(polygon, errorMsg)) {
            cerr << "Ошибка: " << errorMsg << endl;
            return 1;
        }

        vector<Point> triangles = triangulate(polygon);
        printTriangles(triangles);
    } catch (const exception& e) {
        cerr << "Ошибка: " << e.what() << endl;
        return 1;
    }
    return 0;
}
