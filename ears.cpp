#include <iostream>
#include <vector>
#include <cmath>
#include <list>
#include <sstream>
#include <fstream>
#include <SFML/Graphics.hpp>

using namespace std;

const double EPS = 1e-6;

struct Point {
    double x, y;
    Point(double x = 0, double y = 0) : x(x), y(y) {}
    Point operator-(const Point& other) const {
        return Point(x - other.x, y - other.y);
    }
    bool operator==(const Point& other) const {
        return fabs(x - other.x) < EPS && fabs(y - other.y) < EPS;
    }
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
        return (min1x <= max2x + EPS && max1x >= min2x - EPS &&
                min1y <= max2y + EPS && max1y >= min2y - EPS);
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

vector<Point> removeDuplicates(const vector<Point>& input) {
    vector<Point> result;
    for (const Point& p : input) {
        bool duplicate = false;
        for (const Point& existing : result) {
            if (existing == p) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) {
            result.push_back(p);
        }
    }
    return result;
}

bool validatePolygon(vector<Point>& polygon, string& errorMsg) {
    polygon = removeDuplicates(polygon);
    int n = polygon.size();

    if (n < 3) {
        errorMsg = "Многоугольник должен содержать как минимум 3 вершины";
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
        cout << "Направление обхода изменено на против часовой стрелки" << endl;
    }
    return true;
}

struct Vertex {
    Point p;
    int index;
    Vertex* prev;
    Vertex* next;
    bool removed;
    bool isConvex;
    bool isEar;
    list<Vertex*>::iterator earIt;
    list<Vertex*>::iterator concaveIt;

    Vertex(const Point& pt, int idx) : p(pt), index(idx), prev(nullptr), next(nullptr),
                                        removed(false), isConvex(false), isEar(false) {}
};

class Triangulator {
private:
    vector<Vertex*> input_points;
    Vertex* head;
    int remaining;
    list<Vertex*> concave_list;
    list<Vertex*> ear_list;

    bool isConvex(Vertex* v) {
        if (!v || v->removed || !v->prev || !v->next || v->prev->removed || v->next->removed) {
            return false;
        }
        return crossProduct(v->prev->p, v->p, v->next->p) >= -EPS;
    }

    bool isEar(Vertex* v) {
        if (!v || v->removed || !v->isConvex) return false;
        if (!v->prev || !v->next || v->prev->removed || v->next->removed) return false;

        for (Vertex* concave : concave_list) {
            if (!concave || concave->removed) continue;
            if (concave == v || concave == v->prev || concave == v->next) continue;
            if (pointInTriangle(concave->p, v->prev->p, v->p, v->next->p)) {
                return false;
            }
        }
        return true;
    }

    void updateVertex(Vertex* v) {
        if (!v || v->removed) return;

        bool wasConvex = v->isConvex;
        v->isConvex = isConvex(v);

        if (!v->isConvex && wasConvex) {
            concave_list.push_back(v);
            v->concaveIt = --concave_list.end();
        } else if (v->isConvex && !wasConvex) {
            if (v->concaveIt != concave_list.end()) {
                concave_list.erase(v->concaveIt);
                v->concaveIt = concave_list.end();
            }
        }

        bool wasEar = v->isEar;
        v->isEar = isEar(v);

        if (v->isEar && !wasEar) {
            ear_list.push_back(v);
            v->earIt = --ear_list.end();
        } else if (!v->isEar && wasEar) {
            if (v->earIt != ear_list.end()) {
                ear_list.erase(v->earIt);
                v->earIt = ear_list.end();
            }
        }
    }

    void removeVertex(Vertex* v) {
        if (!v || v->removed) return;

        if (v->isEar && v->earIt != ear_list.end()) {
            ear_list.erase(v->earIt);
            v->earIt = ear_list.end();
        }

        if (!v->isConvex && v->concaveIt != concave_list.end()) {
            concave_list.erase(v->concaveIt);
            v->concaveIt = concave_list.end();
        }

        Vertex* prev = v->prev;
        Vertex* next = v->next;

        if (prev) prev->next = next;
        if (next) next->prev = prev;

        v->removed = true;
        v->prev = nullptr;
        v->next = nullptr;

        if (head == v) {
            head = next;
        }

        remaining--;

        if (prev && !prev->removed) updateVertex(prev);
        if (next && !next->removed) updateVertex(next);
    }

    void initializeLists() {
        if (!head || remaining < 3) return;

        Vertex* current = head;
        int count = 0;
        while (count < remaining && current != nullptr) {
            current->isConvex = isConvex(current);

            if (!current->isConvex) {
                concave_list.push_back(current);
                current->concaveIt = --concave_list.end();
            }

            current->isEar = isEar(current);

            if (current->isEar) {
                ear_list.push_back(current);
                current->earIt = --ear_list.end();
            }

            current = current->next;
            count++;
            if (current == head || current == nullptr) break;
        }
    }

public:
    Triangulator(const vector<Point>& polygon) {
        int n = polygon.size();
        input_points.resize(n);
        for (int i = 0; i < n; ++i) {
            input_points[i] = new Vertex(polygon[i], i);
        }
        for (int i = 0; i < n; ++i) {
            input_points[i]->prev = input_points[(i - 1 + n) % n];
            input_points[i]->next = input_points[(i + 1) % n];
            input_points[i]->concaveIt = concave_list.end();
            input_points[i]->earIt = ear_list.end();
        }
        head = input_points[0];
        remaining = n;
    }

    ~Triangulator() {
        for (Vertex* v : input_points) {
            delete v;
        }
    }

    vector<Point> triangulate() {
        vector<Point> triangles;
        if (remaining < 3) return triangles;

        initializeLists();

        while (remaining > 3) {
            if (ear_list.empty()) {
                cerr << "Ошибка: список ушей пуст" << endl;
                break;
            }

            Vertex* v_curr = nullptr;
            list<Vertex*>::iterator it;

            for (it = ear_list.begin(); it != ear_list.end(); ++it) {
                Vertex* v = *it;
                if (v && !v->removed && v->isEar) {
                    v_curr = v;
                    break;
                }
            }

            if (!v_curr) {
                cerr << "Ошибка: не найдено валидное ухо" << endl;
                break;
            }

            ear_list.erase(it);
            v_curr->earIt = ear_list.end();

            Vertex* v_prev = v_curr->prev;
            Vertex* v_next = v_curr->next;

            if (!v_prev || !v_next || v_prev->removed || v_next->removed) {
                continue;
            }

            triangles.push_back(v_prev->p);
            triangles.push_back(v_curr->p);
            triangles.push_back(v_next->p);

            removeVertex(v_curr);
        }

        if (remaining == 3) {
            vector<Vertex*> remaining_vertices;
            Vertex* current = head;
            int count = 0;
            while (count < (int)input_points.size() && remaining_vertices.size() < 3) {
                if (current && !current->removed) {
                    remaining_vertices.push_back(current);
                }
                current = current ? current->next : nullptr;
                count++;
                if (current == head) break;
            }

            if (remaining_vertices.size() == 3) {
                triangles.push_back(remaining_vertices[0]->p);
                triangles.push_back(remaining_vertices[1]->p);
                triangles.push_back(remaining_vertices[2]->p);
            }
        }

        return triangles;
    }
};

vector<Point> readPointsFromArgs(int argc, char* argv[]) {
    vector<Point> points;

    if (argc == 2) {
        ifstream file(argv[1]);
        if (!file.is_open()) {
            throw runtime_error("Не удалось открыть файл: " + string(argv[1]));
        }

        string content;
        string line;
        while (getline(file, line)) {
            content += line;
        }
        file.close();

        stringstream ss(content);
        string token;
        vector<double> coords;

        while (getline(ss, token, ',')) {
            if (!token.empty()) {
                double val;
                stringstream valStream(token);
                if (valStream >> val) {
                    coords.push_back(val);
                }
            }
        }

        for (size_t i = 0; i + 1 < coords.size(); i += 2) {
            points.push_back(Point(coords[i], coords[i + 1]));
        }
    } else {
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
    }

    return points;
}

void visualize(const vector<Point>& polygon, const vector<Point>& triangles) {
    if (polygon.empty() || triangles.empty()) {
        cerr << "Нет данных для визуализации" << endl;
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

    double range = max(maxX - minX, maxY - minY);
    double centerX = (minX + maxX) / 2;
    double centerY = (minY + maxY) / 2;

    const int WINDOW_WIDTH = 800;
    const int WINDOW_HEIGHT = 800;
    const int MARGIN = 40;

    auto toScreenX = [&](double x) {
        return (x - centerX) / range * (WINDOW_WIDTH - 2 * MARGIN) / 2 + WINDOW_WIDTH / 2;
    };
    auto toScreenY = [&](double y) {
        return -(y - centerY) / range * (WINDOW_HEIGHT - 2 * MARGIN) / 2 + WINDOW_HEIGHT / 2;
    };

    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Triangulation Visualization");

    sf::Font font;
    std::string fontPaths[] = {
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf"
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

        for (size_t i = 0; i + 2 < triangles.size(); i += 3) {
            sf::ConvexShape triangle;
            triangle.setPointCount(3);
            triangle.setPoint(0, sf::Vector2f(toScreenX(triangles[i].x), toScreenY(triangles[i].y)));
            triangle.setPoint(1, sf::Vector2f(toScreenX(triangles[i+1].x), toScreenY(triangles[i+1].y)));
            triangle.setPoint(2, sf::Vector2f(toScreenX(triangles[i+2].x), toScreenY(triangles[i+2].y)));
            triangle.setFillColor(sf::Color(200, 220, 255, 180));
            triangle.setOutlineThickness(0);
            window.draw(triangle);
        }

        for (size_t i = 0; i < polygon.size(); ++i) {
            size_t j = (i + 1) % polygon.size();
            sf::Vertex line[] = {
                sf::Vertex(sf::Vector2f(toScreenX(polygon[i].x), toScreenY(polygon[i].y)), sf::Color(255, 0, 0)),
                sf::Vertex(sf::Vector2f(toScreenX(polygon[j].x), toScreenY(polygon[j].y)), sf::Color(255, 0, 0))
            };
            line[0].color = sf::Color(255, 0, 0, 255);
            line[1].color = sf::Color(255, 0, 0, 255);
            window.draw(line, 2, sf::Lines);
        }

        for (size_t i = 0; i + 2 < triangles.size(); i += 3) {
            for (int j = 0; j < 3; ++j) {
                int k = (j + 1) % 3;
                Point p1 = triangles[i + j];
                Point p2 = triangles[i + k];
                bool isBorder = false;
                for (size_t v = 0; v < polygon.size(); ++v) {
                    size_t next = (v + 1) % polygon.size();
                    if ((p1 == polygon[v] && p2 == polygon[next]) ||
                        (p1 == polygon[next] && p2 == polygon[v])) {
                        isBorder = true;
                        break;
                    }
                }
                if (!isBorder) {
                    sf::Vertex line[] = {
                        sf::Vertex(sf::Vector2f(toScreenX(p1.x), toScreenY(p1.y)), sf::Color(100, 150, 255, 120)),
                        sf::Vertex(sf::Vector2f(toScreenX(p2.x), toScreenY(p2.y)), sf::Color(100, 150, 255, 120))
                    };
                    window.draw(line, 2, sf::Lines);
                }
            }
        }
        int step = 1;
        if (polygon.size() > 200) step = polygon.size() / 100 + 1;

        for (size_t i = 0; i < polygon.size(); i += step) {
            const Point& p = polygon[i];
            sf::CircleShape circle(3);
            circle.setPosition(toScreenX(p.x) - 3, toScreenY(p.y) - 3);
            circle.setFillColor(sf::Color(255, 0, 0));
            window.draw(circle);

            if (fontLoaded && i % (step * 5) == 0 && i < polygon.size()) {
                sf::Text text;
                text.setFont(font);
                text.setString(to_string(i));
                text.setCharacterSize(10);
                text.setFillColor(sf::Color::Black);
                text.setPosition(toScreenX(p.x) + 5, toScreenY(p.y) - 5);
                window.draw(text);
            }
        }

        if (fontLoaded) {
            sf::Text info;
            info.setFont(font);
            info.setString("Vertices: " + to_string(polygon.size()) +
                          "  Triangles: " + to_string(triangles.size() / 3) +
                          "\nESC - exit\nRed: boundary | Blue: diagonals");
            info.setCharacterSize(14);
            info.setFillColor(sf::Color::Black);
            info.setPosition(10, 10);
            window.draw(info);
        }

        window.display();
    }
}

int main(int argc, char* argv[]) {
    try {
        vector<Point> polygon;
        string errorMsg;

        if (argc < 2) {
            cerr << "Использование: " << argv[0] << " [файл.txt] или x1 y1 x2 y2 ..." << endl;
            return 1;
        }

        polygon = readPointsFromArgs(argc, argv);

        cout << "Входной многоугольник (" << polygon.size() << " вершин)" << endl;

        if (!validatePolygon(polygon, errorMsg)) {
            cerr << "Ошибка: " << errorMsg << endl;
            return 1;
        }

        cout << "После удаления дубликатов (" << polygon.size() << " вершин)" << endl;

        Triangulator triangulator(polygon);
        vector<Point> triangles = triangulator.triangulate();

        cout << "Треугольников: " << triangles.size() / 3 << endl;

        if (!triangles.empty()) {
            visualize(polygon, triangles);
        } else {
            cerr << "Предупреждение: триангуляция не выполнена" << endl;
        }

    } catch (const exception& e) {
        cerr << "Ошибка: " << e.what() << endl;
        return 1;
    }
    return 0;
}
