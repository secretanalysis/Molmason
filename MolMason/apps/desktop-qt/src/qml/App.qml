import QtQuick 2.15
import QtQuick.Controls 2.15
import QtQuick.Layouts 1.15

ApplicationWindow {
    id: root
    visible: true
    width: 1200
    height: 800
    title: "MolMason"
    
    RowLayout {
        anchors.fill: parent
        spacing: 0
        
        // Sidebar
        Rectangle {
            Layout.preferredWidth: 200
            Layout.fillHeight: true
            color: "#2d2d2d"
            
            ColumnLayout {
                anchors.fill: parent
                anchors.margins: 10
                
                Label {
                    text: "MolMason"
                    font.pixelSize: 24
                    font.bold: true
                    color: "white"
                }
                
                Label {
                    text: "Local-first molecule design"
                    font.pixelSize: 12
                    color: "#aaa"
                }
                
                Item { Layout.preferredHeight: 20 }
                
                Repeater {
                    model: ["Explorer", "Protocol", "Feedback", "Diagnostics"]
                    
                    Button {
                        Layout.fillWidth: true
                        text: modelData
                        flat: true
                        contentItem: Text {
                            text: parent.text
                            color: "white"
                            horizontalAlignment: Text.AlignLeft
                        }
                    }
                }
                
                Item { Layout.fillHeight: true }
            }
        }
        
        // Main content
        Rectangle {
            Layout.fillWidth: true
            Layout.fillHeight: true
            color: "#1e1e1e"
            
            Label {
                anchors.centerIn: parent
                text: "MolMason v2.0.0"
                font.pixelSize: 32
                color: "#666"
            }
        }
    }
}
