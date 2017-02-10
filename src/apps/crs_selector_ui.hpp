#ifndef __CRS_SELECTOR_UI_HPP__
#define __CRS_SELECTOR_UI_HPP__

#include <map>
#include <string>

#include <QWidget>
#include <QDir>
#include <QCompleter>
#include <QString>
#include <QKeyEvent>
#include <QDialog>

#include "ui_crs_selector.h"

namespace geo {

	namespace ui {

		class CRSSelector: public QDialog, public Ui::CRSSelector {
			Q_OBJECT
		private:
			std::map<int, std::string> m_hcrs;
			std::map<int, std::string> m_vcrs;
			QCompleter *m_vcomp;
			QCompleter *m_hcomp;
			QStringList m_vlst;
			QStringList m_hlst;
			int m_vsrid;
			int m_hsrid;
			void initUi();
			void loadCrs(std::map<int, std::string> &target,
					const std::string &filename);
			void updateFields();
		protected:
			void keyPressEvent(QKeyEvent *e);
		public:
			CRSSelector(QWidget *parent = Q_NULLPTR, Qt::WindowFlags f =
					Qt::WindowFlags());
			void enableVertical(bool);
			void enableHorizontal(bool);
			int getHorizontalSRID();
			int getVerticalSRID();
			void setHorizontalSRID(int);
			void setVerticalSRID(int);
		public slots:
			void textUpdate();
			void cancelClicked();
			void selectClicked();
		};

	}

}

#endif

