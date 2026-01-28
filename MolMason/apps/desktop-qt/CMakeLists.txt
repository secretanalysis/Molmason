cmake_minimum_required(VERSION 3.16)
project(MolMason VERSION 2.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Qt6 REQUIRED COMPONENTS Core Quick Widgets)

qt_add_executable(molmason
    src/main.cpp
)

qt_add_qml_module(molmason
    URI MolMason
    VERSION 1.0
    QML_FILES
        src/qml/App.qml
        src/qml/screens/Explorer.qml
        src/qml/screens/Protocol.qml
        src/qml/screens/Feedback.qml
        src/qml/screens/Diagnostics.qml
)

target_link_libraries(molmason PRIVATE
    Qt6::Core
    Qt6::Quick
    Qt6::Widgets
)
