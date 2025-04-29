from .main_window import QApplication, MainWindow
import sys 

app = QApplication(sys.argv)
w = MainWindow(); w.show()
sys.exit(app.exec()) 
