#include <string>
#include <vector>
#include <map>

#include <QtWidgets/QApplication>
#include <QtWidgets/QMessageBox>

#include <cpl_conv.h>

#include "geo.hpp"
#include "treetops.hpp"
#include "ui/treetops_ui.hpp"

#pragma comment(linker, "/SUBSYSTEM:windows /ENTRY:mainCRTStartup")

int runWithGui(int argc, char** argv) {

	class tt_app : public QApplication {
	public:
		
	tt_app(int& argc, char** argv) : QApplication(argc, argv) {}
		
		bool notify(QObject* receiver, QEvent* e) {
			try {
				return QApplication::notify(receiver, e);
			} catch(const std::exception &ex) {
				QMessageBox err;
				err.setText("Error");
				err.setInformativeText(QString(ex.what()));
				err.exec();
				return false;
			} catch (...) {
				g_warn("Some qt exception.")
			}
			return false;
		}
	};

	// These must be set before TTForm is constructed because
	// QSettings needs them for context.
	QCoreApplication::setOrganizationName("Rob Skelly");
	QCoreApplication::setOrganizationDomain("dijital.ca");
	QCoreApplication::setApplicationName("Treetops");

	tt_app q(argc, argv);
	geo::ui::TTForm f;
	f.showForm();
	return q.exec();
}

int main(int argc, char** argv) {

	try {

		return runWithGui(argc, argv);

	} catch (const std::exception &e) {
		std::cerr << e.what() << std::endl;
		return 1;
	}

}

