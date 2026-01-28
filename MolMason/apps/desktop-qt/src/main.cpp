#include <QGuiApplication>
#include <QQmlApplicationEngine>

int main(int argc, char *argv[]) {
    QGuiApplication app(argc, argv);
    app.setApplicationName("MolMason");
    app.setApplicationVersion("2.0.0");
    app.setOrganizationName("KingCermak");
    
    QQmlApplicationEngine engine;
    engine.load(QUrl(QStringLiteral("qrc:/MolMason/src/qml/App.qml")));
    
    return app.exec();
}
