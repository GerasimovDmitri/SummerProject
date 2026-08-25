#include <iostream>
#include <vector>
#include <cmath>
#include <list>
#include <sstream>
#include <SFML/Graphics.hpp>

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
    int index;
    Vertex* prev;
    Vertex* next;
    bool removed;
    bool isConvex;
    bool isEar;
    list<Vertex*>::iterator earIt;

    Vertex(const Point& pt, int idx) : p(pt), index(idx), prev(nullptr), next(nullptr),
                                        removed(false), isConvex(false), isEar(false) {}
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

bool pointInTriangle(const Point& p, const Point& a, const Point& b, const Point& c) {
    double d1 = crossProduct(a, b, p);
    double d2 = crossProduct(b, c, p);
    double d3 = crossProduct(c, a, p);
    bool hasNeg = (d1 < -EPS) || (d2 < -EPS) || (d3 < -EPS);
    bool hasPos = (d1 > EPS) || (d2 > EPS) || (d3 > EPS);
    return !(hasNeg && hasPos);
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

class Triangulator {
private:
    vector<Vertex*> vertices;
    Vertex* head;
    int remaining;
    list<Vertex*> earList;

    bool isConvex(Vertex* v) {
        return crossProduct(v->prev->p, v->p, v->next->p) >= -EPS;
    }

    bool isEar(Vertex* v) {
        if (!v->isConvex) return false;
        Vertex* current = v->next->next;
        while (current != v->prev) {
            if (!current->removed && pointInTriangle(current->p, v->prev->p, v->p, v->next->p)) {
                return false;
            }
            current = current->next;
        }
        return true;
    }

    void updateVertex(Vertex* v) {
        if (v->removed) return;
        v->isConvex = isConvex(v);
        v->isEar = isEar(v);
    }

    void addToEarList(Vertex* v) {
        if (v->removed || !v->isEar) return;
        earList.push_back(v);
        v->earIt = --earList.end();
    }

    void removeFromEarList(Vertex* v) {
        if (!v->isEar) return;
        earList.erase(v->earIt);
        v->isEar = false;
    }

    void removeVertex(Vertex* v) {
        Vertex* prev = v->prev;
        Vertex* next = v->next;

        prev->next = next;
        next->prev = prev;

        v->removed = true;
        removeFromEarList(v);

        if (head == v) {
            head = next;
        }

        remaining--;

        removeFromEarList(prev);
        removeFromEarList(next);
        
        updateVertex(prev);
        updateVertex(next);
        
        addToEarList(prev);
        addToEarList(next);
    }

    void initializeEarList() {
        Vertex* current = head;
        for (int i = 0; i < remaining; ++i) {
            current->isConvex = isConvex(current);
            current->isEar = isEar(current);
            if (current->isEar) {
                earList.push_back(current);
                current->earIt = --earList.end();
            }
            current = current->next;
        }
    }

public:
    Triangulator(const vector<Point>& polygon) {
        int n = polygon.size();
        vertices.resize(n);
        for (int i = 0; i < n; ++i) {
            vertices[i] = new Vertex(polygon[i], i);
        }
        for (int i = 0; i < n; ++i) {
            vertices[i]->prev = vertices[(i - 1 + n) % n];
            vertices[i]->next = vertices[(i + 1) % n];
        }
        head = vertices[0];
        remaining = n;
    }

    ~Triangulator() {
        for (Vertex* v : vertices) {
            delete v;
        }
    }

    vector<Point> triangulate() {
        vector<Point> triangles;
        if (remaining < 3) return triangles;

        initializeEarList();

        while (remaining > 3) {
            if (earList.empty()) {
                cerr << "Предупреждение: ухо не найдено, завершение" << endl;
                break;
            }

            Vertex* ear = earList.front();
            earList.pop_front();

            if (ear->removed || !ear->isEar) {
                continue;
            }

            Vertex* prev = ear->prev;
            Vertex* next = ear->next;

            triangles.push_back(prev->p);
            triangles.push_back(ear->p);
            triangles.push_back(next->p);

            removeVertex(ear);
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

        return triangles;
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

void visualize(const vector<Point>& polygon, const vector<Point>& triangles) {
    if (triangles.empty()) {
        cerr << "Триангуляция пуста, визуализация невозможна" << endl;
        return;
    }

    double minX = polygon[0].x, maxX = polygon[0].x;
    double minY = polygon[0].y, maxY = polygon[0].y;
    for (const Point& p : polygon) {
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
    const int WINDOW_HEIGHT = 800;

    auto toScreenX = [&](double x) {
        return (x - minX) / (maxX - minX) * (WINDOW_WIDTH - 40) + 20;
    };
    auto toScreenY = [&](double y) {
        return (maxY - y) / (maxY - minY) * (WINDOW_HEIGHT - 40) + 20;
    };

    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Triangulation Visualization");

    sf::Font font;
    std::string fontPaths[] = {
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
        "/usr/share/fonts/truetype/arial/arial.ttf",
        "/System/Library/Fonts/Helvetica.ttf",
        "/Windows/Fonts/arial.ttf"
    };

    bool fontLoaded = false;
    for (const auto& path : fontPaths) {
        if (font.loadFromFile(path)) {
            fontLoaded = true;
            break;
        }
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

        for (size_t i = 0; i < triangles.size(); i += 3) {
            Point p1 = triangles[i];
            Point p2 = triangles[i + 1];
            Point p3 = triangles[i + 2];

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

        for (size_t i = 0; i < polygon.size(); ++i) {
            size_t j = (i + 1) % polygon.size();
            sf::Vertex line[] = {
                sf::Vertex(sf::Vector2f(toScreenX(polygon[i].x), toScreenY(polygon[i].y)), sf::Color::Red),
                sf::Vertex(sf::Vector2f(toScreenX(polygon[j].x), toScreenY(polygon[j].y)), sf::Color::Red)
            };
            window.draw(line, 2, sf::Lines);
        }

        for (size_t i = 0; i < polygon.size(); ++i) {
            const Point& p = polygon[i];
            sf::CircleShape circle(6);
            circle.setPosition(toScreenX(p.x) - 6, toScreenY(p.y) - 6);
            circle.setFillColor(sf::Color::Red);
            circle.setOutlineColor(sf::Color::Black);
            circle.setOutlineThickness(1);
            window.draw(circle);

            if (fontLoaded) {
                sf::Text text;
                text.setFont(font);
                text.setString(to_string(i));
                text.setCharacterSize(14);
                text.setFillColor(sf::Color::Black);
                text.setPosition(toScreenX(p.x) + 8, toScreenY(p.y) - 8);
                window.draw(text);
            }
        }

        if (fontLoaded) {
            sf::Text info;
            info.setFont(font);
            info.setString("Vertices: " + to_string(polygon.size()) +
                          "  Triangles: " + to_string(triangles.size() / 3) +
                          "\nESC - exit");
            info.setCharacterSize(16);
            info.setFillColor(sf::Color::Black);
            info.setPosition(10, 10);
            window.draw(info);
        }

        window.display();
    }
}

int main(int argc, char* argv[]) {
    try {
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

        Triangulator triangulator(polygon);
        vector<Point> triangles = triangulator.triangulate();
        printTriangles(triangles);
        visualize(polygon, triangles);

    } catch (const exception& e) {
        cerr << "Ошибка: " << e.what() << endl;
        return 1;
    }
    return 0;
}
