#ifndef __LASGRID_UI_HPP__
#define __LASGRID_UI_HPP__

#include <set>

#include <QtWidgets/QWidget>
#include <QDir>

#include "ui_lasgrid.h"

namespace geotools {

	namespace ui {

		class MosaicForm : public QWidget, public Ui::LasgridForm {
			Q_OBJECT
		private:
			QWidget *m_form;
			std::string m_destFile;
			std::vector<std::string> m_tifFiles;
			QDir *m_last;
			int m_threads;
			double m_distance;
			int m_tileSize;
			void updateFileList();
			void updateFileButtons();
			void checkRun();
		public:
			MosaicForm(QWidget *p = Q_NULLPTR);
			void setupUi(QWidget *Form);
			~MosaicForm();
			void setFileStatus(float);
			void setOverallStatus(float);

		public slots:
			void fileListSelectionChanged();
			void selectFilesClicked();
			void removeFilesClicked();
			void clearFilesClicked();
			void cancelClicked();
			void runClicked();
			void destFileClicked();
			void upClicked();
			void downClicked();
			void distanceChanged(double);
			void threadsChanged(int);
			void tileSizeChanged(int);
		};

	}

}

#endif

