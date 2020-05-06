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

#include <QtCore/QDir>
#include <QtCore/QObject>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QWidget>

namespace geo {

namespace ui {

class FileList: public QObject {
	Q_OBJECT
private:
	std::vector<std::string> m_files;			///<! The internal file list.
	QPushButton *m_btnAddFiles;					///<! The add files button.
	QPushButton *m_btnRemoveSelectedFiles;		///<! The remove selected files button.
	QPushButton *m_btnRemoveAllFiles;			///<! The remove all files button.
	QListWidget *m_lstFiles;					///<! The file list widget.
	QWidget *m_parent;							///<! The parent widget.
	QDir m_last;								///<! The last-used directory.
	QString m_filter;							///<! The file name filter.

	/// Updates the file list state.
	void updateFileList();
	/// Updates the button states.
	void updateButtons();

public:
	/// Construct an empty file list.
	FileList();

	/**
	 * Initialize the file list.
	 *
	 * \param parent The parent widget.
	 * \param btnAddFiles The add files button.
	 * \param btnRemoveAllFiles The remove all files button.
	 * \param btnRemoveSelectedFiles The remove selected files button.
	 * \param lstFiles the list widget containing the files.
	 * \param last The last-used directory.
	 * \param filter A string filter to use on the file list.
	 */
	void init(QWidget* parent, QPushButton* btnAddFiles,
			QPushButton* btnRemoveAllFiles, QPushButton* btnRemoveSelectedFiles,
			QListWidget* lstFiles, QDir& last, QString& filter);

	/**
	 * Return a vector containing the files in the list.
	 *
	 * \return A vector containing the file paths in the list.
	 */
	std::vector<std::string> files();

	/**
	 * Populate the file list with the given files.
	 *
	 * \param files A vector containing file paths.
	 */
	void setFiles(const std::vector<std::string>& files);

	virtual ~FileList();

public slots:
	/// Called when the selected file is changed.
	void fileListSelectionChanged();
	/// Called with the file removal button is clicked.
	void removeSelectedFilesClicked();
	/// Called with the remove all button is clicked.
	void removeAllFilesClicked();
	/// Called when the add file button is clicked.
	void addFilesClicked();

signals:
	/// Signals that the file list has changed.
	void fileListChanged();

};

}
}

#endif /* FILELIST_HPP */

