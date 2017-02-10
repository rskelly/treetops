/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   filelist.hpp
 * Author: rob
 *
 * Created on November 27, 2016, 1:04 PM
 */

#ifndef FILELIST_HPP
#define FILELIST_HPP

#include <QObject>
#include <QPushButton>
#include <QListWidget>
#include <QWidget>
#include <QDir>

namespace geo {

namespace ui {

class FileList: public QObject {
	Q_OBJECT
private:
	std::vector<std::string> m_files;
	QPushButton *m_btnAddFiles;
	QPushButton *m_btnRemoveSelectedFiles;
	QPushButton *m_btnRemoveAllFiles;
	QListWidget *m_lstFiles;
	QWidget *m_parent;
	QDir m_last;
	QString m_filter;

	void updateFileList();
	void updateButtons();

public:
	FileList();

	void init(QWidget *parent, QPushButton *btnAddFiles,
			QPushButton *btnRemoveAllFiles, QPushButton *btnRemoveSelectedFiles,
			QListWidget *lstFiles, QDir &last, QString &filter);

	std::vector<std::string> files();
	void setFiles(const std::vector<std::string> &files);

	virtual ~FileList();

public slots:
	void fileListSelectionChanged();
	void removeSelectedFilesClicked();
	void removeAllFilesClicked();
	void addFilesClicked();

	signals:
	void fileListChanged();

};

}
}

#endif /* FILELIST_HPP */

