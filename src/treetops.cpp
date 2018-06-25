#include <string>
#include <vector>
#include <map>

#include <QtWidgets/QApplication>
#include <QtWidgets/QMessageBox>

#include "geo.hpp"
#include "treetops.hpp"

#ifdef WITH_GUI
#include "treetops_ui.hpp"
#include "cpl_conv.h"
#endif

#pragma comment(linker, "/SUBSYSTEM:windows /ENTRY:mainCRTStartup")

int runWithGui(int argc, char **argv) {
#ifdef WITH_GUI

#ifdef _MSC_VER
	CPLSetConfigOption("GDAL_DATA", "./gdal-data");
#endif

	class TTApplication : public QApplication {
	public:
		TTApplication(int &argc, char **argv) : QApplication(argc, argv) {}
		bool notify(QObject *receiver, QEvent *e) {
			try {
				return QApplication::notify(receiver, e);
			} catch(const std::exception &ex) {
				QMessageBox err;
				err.setText("Error");
				err.setInformativeText(QString(ex.what()));
				err.exec();
				return false;
			}
		}
	};

	QCoreApplication::setOrganizationName("Rob Skelly");
	QCoreApplication::setOrganizationDomain("dijital.ca");
	QCoreApplication::setApplicationName("Treetops");

	TTApplication q(argc, argv);
	geo::ui::TreetopsForm f;
	f.showForm();
	return q.exec();
#else
	std::cerr << "GUI not enabled." << std::endl;
	return 1;
#endif
}

int main(int argc, char **argv) {

	using namespace geo::treetops;
	using namespace geo::treetops::util;
	using namespace geo::treetops::config;

	try {

		return runWithGui(argc, argv);

	} catch (const std::exception &e) {
		std::cerr << e.what() << std::endl;
		return 1;
	}

	return 0;
}

