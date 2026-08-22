#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>
#include <sstream>
#include <SFML/Graphics.hpp>

using namespace std;

const double EPS = 1e-6;

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
    
    bool isOnBoundary(const Point& p) const {
        double dist = sqrt(center.dist2(p));
        return fabs(dist - radius) < 0.01;
    }
};

class MinimumEnclosingCircle {
private:
    vector<Point> points;
    vector<int> indices;

    Circle circleFrom2Points(const Point& p1, const Point& p2) {
        Point center((p1.x + p2.x) / 2, (p1.y + p2.y) / 2);
        double radius = sqrt(center.dist2(p1));
        return Circle(center, radius);
    }

    Circle circleFrom3Points(const Point& p1, const Point& p2, const Point& p3) {
        double d = 2 * (p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y));
        if (fabs(d) < EPS) {
            return circleFrom2Points(p1, p2);
        }
        double ux = ((p1.x*p1.x + p1.y*p1.y) * (p2.y - p3.y) +
                     (p2.x*p2.x + p2.y*p2.y) * (p3.y - p1.y) +
                     (p3.x*p3.x + p3.y*p3.y) * (p1.y - p2.y)) / d;
        double uy = -((p1.x*p1.x + p1.y*p1.y) * (p2.x - p3.x) +
                      (p2.x*p2.x + p2.y*p2.y) * (p3.x - p1.x) +
                      (p3.x*p3.x + p3.y*p3.y) * (p1.x - p2.x)) / d;
        Point center(ux, uy);
        double radius = sqrt(center.dist2(p1));
        return Circle(center, radius);
    }

public:
    MinimumEnclosingCircle(const vector<Point>& pts) : points(pts) {
        indices.resize(points.size());
        for (size_t i = 0; i < points.size(); i++) indices[i] = i;
    }

    Circle compute() {
        if (points.empty()) return Circle(Point(0, 0), 0);
        if (points.size() == 1) return Circle(points[0], 0);
        
        unsigned seed = chrono::steady_clock::now().time_since_epoch().count();
        shuffle(indices.begin(), indices.end(), mt19937(seed));

        Circle c(points[indices[0]], 0);
        
        for (size_t i = 1; i < points.size(); i++) {
            const Point& p = points[indices[i]];
            if (!c.contains(p)) {
                c = Circle(p, 0);
                for (size_t j = 0; j < i; j++) {
                    const Point& q = points[indices[j]];
                    if (!c.contains(q)) {
                        c = circleFrom2Points(p, q);
                        for (size_t k = 0; k < j; k++) {
                            const Point& r = points[indices[k]];
                            if (!c.contains(r)) {
                                c = circleFrom3Points(p, q, r);
                            }
                        }
                    }
                }
            }
        }
        
        return c;
    }
    
    int countPointsOnBoundary(const Circle& c) const {
        int count = 0;
        for (size_t i = 0; i < points.size(); i++) {
            if (c.isOnBoundary(points[i])) {
                count++;
            }
        }
        return count;
    }
};

void visualize(const vector<Point>& points, const Circle& minCircle) {
    if (points.empty()) {
        cerr << "No points to visualize" << endl;
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

    double margin = max(maxX - minX, maxY - minY) * 0.3 + minCircle.radius * 1.2 + 2.0;
    if (margin < 5.0) margin = 5.0;
    
    minX -= margin;
    maxX += margin;
    minY -= margin;
    maxY += margin;

    const int WINDOW_SIZE = 800;
    const int WINDOW_WIDTH = WINDOW_SIZE;
    const int WINDOW_HEIGHT = WINDOW_SIZE;

    double rangeX = maxX - minX;
    double rangeY = maxY - minY;
    double maxRange = max(rangeX, rangeY);
    
    double centerX = (minX + maxX) / 2.0;
    double centerY = (minY + maxY) / 2.0;
    
    double newMinX = centerX - maxRange / 2.0;
    double newMaxX = centerX + maxRange / 2.0;
    double newMinY = centerY - maxRange / 2.0;
    double newMaxY = centerY + maxRange / 2.0;

    auto toScreenX = [&](double x) {
        return (x - newMinX) / (newMaxX - newMinX) * (WINDOW_WIDTH - 40) + 20;
    };
    auto toScreenY = [&](double y) {
        return (newMaxY - y) / (newMaxY - newMinY) * (WINDOW_HEIGHT - 40) + 20;
    };

    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Minimum Enclosing Circle");
    
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

        if (minCircle.radius > EPS) {
            double screenRadius = minCircle.radius * (WINDOW_WIDTH - 40) / (newMaxX - newMinX);
            
            sf::CircleShape circle;
            circle.setRadius(screenRadius);
            circle.setOrigin(screenRadius, screenRadius);
            circle.setPosition(toScreenX(minCircle.center.x), toScreenY(minCircle.center.y));
            circle.setFillColor(sf::Color::Transparent);
            circle.setOutlineColor(sf::Color(255, 0, 0, 255));
            circle.setOutlineThickness(3.0f);
            window.draw(circle);
            
            sf::CircleShape centerDot(6);
            centerDot.setOrigin(6, 6);
            centerDot.setPosition(toScreenX(minCircle.center.x), toScreenY(minCircle.center.y));
            centerDot.setFillColor(sf::Color(255, 0, 0));
            centerDot.setOutlineColor(sf::Color::Black);
            centerDot.setOutlineThickness(1.5);
            window.draw(centerDot);
        }

        for (size_t i = 0; i < points.size(); i++) {
            const Point& p = points[i];
            
            double dist = sqrt(minCircle.center.dist2(p));
            bool onBoundary = fabs(dist - minCircle.radius) < 0.01;
            
            sf::CircleShape circle(8);
            circle.setOrigin(8, 8);
            circle.setPosition(toScreenX(p.x), toScreenY(p.y));
            
            if (onBoundary) {
                circle.setFillColor(sf::Color(255, 0, 0));
                circle.setOutlineColor(sf::Color::Black);
                circle.setOutlineThickness(2.5);
            } else {
                circle.setFillColor(sf::Color(0, 0, 255));
                circle.setOutlineColor(sf::Color::Black);
                circle.setOutlineThickness(2.0);
            }
            window.draw(circle);

            if (fontLoaded) {
                sf::Text text;
                text.setFont(font);
                text.setString(to_string(i));
                text.setCharacterSize(16);
                text.setFillColor(sf::Color::Black);
                text.setStyle(sf::Text::Bold);
                text.setPosition(toScreenX(p.x) + 14, toScreenY(p.y) - 10);
                window.draw(text);
            }
        }
        
        if (fontLoaded) {
            int boundaryCount = 0;
            for (size_t i = 0; i < points.size(); i++) {
                double dist = sqrt(minCircle.center.dist2(points[i]));
                if (fabs(dist - minCircle.radius) < 0.01) boundaryCount++;
            }
            
            string info = "Points: " + to_string(points.size()) + 
                         "\nCenter: (" + to_string(minCircle.center.x) + ", " + 
                         to_string(minCircle.center.y) + ")" +
                         "\nRadius: " + to_string(minCircle.radius) +
                         "\nPoints on boundary: " + to_string(boundaryCount) +
                         "\nRed points - on boundary" +
                         "\nBlue points - inside" +
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
        }

        if (points.size() < 1) {
            cerr << "Ошибка: нужно минимум 1 точка" << endl;
            cerr << "Использование: " << argv[0] << " x1 y1 x2 y2 x3 y3 ..." << endl;
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
        cout << "  Радиус: " << minCircle.radius << endl;
        cout << "  Точек на границе: " << mec.countPointsOnBoundary(minCircle) << endl << endl;
        visualize(points, minCircle);
        
    } catch (const exception& e) {
        cerr << "Ошибка: " << e.what() << endl;
        return 1;
    }
    return 0;
}
