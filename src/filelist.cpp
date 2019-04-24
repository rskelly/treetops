/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   filelist.cpp
 * Author: rob
 * 
 * Created on November 27, 2016, 1:04 PM
 */

#include <set>

#include <QtCore/QDir>
#include <QtCore/QString>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QWidget>

#include "filelist.hpp"

using namespace geo::ui;

FileList::FileList() :
	m_btnAddFiles(nullptr),
	m_btnRemoveSelectedFiles(nullptr),
	m_btnRemoveAllFiles(nullptr),
	m_lstFiles(nullptr),
	m_parent(nullptr) {
}

void FileList::init(QWidget *parent, QPushButton *btnAddFiles,
		QPushButton *btnRemoveAllFiles, QPushButton *btnRemoveSelectedFiles,
		QListWidget *lstFiles, QDir &last, QString &filter) {
	m_parent = parent;
	m_btnAddFiles = btnAddFiles;
	m_btnRemoveAllFiles = btnRemoveAllFiles;
	m_btnRemoveSelectedFiles = btnRemoveSelectedFiles;
	m_lstFiles = lstFiles;
	m_last = last;
	m_filter = filter;
	connect(btnAddFiles, SIGNAL(clicked()), this, SLOT(addFilesClicked()));
	connect(btnRemoveAllFiles, SIGNAL(clicked()), this,
			SLOT(removeAllFilesClicked()));
	connect(btnRemoveSelectedFiles, SIGNAL(clicked()), this,
			SLOT(removeSelectedFilesClicked()));
	connect(lstFiles, SIGNAL(itemSelectionChanged()), this,
			SLOT(fileListSelectionChanged()));
}

void FileList::updateFileList() {
	while (m_lstFiles->count())
		m_lstFiles->takeItem(0);
	for (const std::string &file : m_files)
		m_lstFiles->addItem(QString(file.c_str()));
	updateButtons();
	emit fileListChanged();
}

void FileList::updateButtons() {
	m_btnRemoveAllFiles->setEnabled(m_lstFiles->count() > 0);
	m_btnRemoveSelectedFiles->setEnabled(
			m_lstFiles->selectedItems().size() > 0);
}

void FileList::fileListSelectionChanged() {
	updateButtons();
}

void FileList::removeSelectedFilesClicked() {
	std::vector<std::string> lst;
	size_t i = 0;
	for (const std::string &file : m_files) {
		QListWidgetItem *item = m_lstFiles->item(i);
		if (!item->isSelected())
			lst.push_back(file);
		++i;
	}
	m_files.clear();
	m_files.assign(lst.begin(), lst.end());
	updateFileList();
}

void FileList::removeAllFilesClicked() {
	m_files.clear();
	updateFileList();
}

void FileList::addFilesClicked() {
	QFileDialog d(m_parent);
	d.setDirectory(m_last);
	d.setFileMode(QFileDialog::ExistingFiles);
	d.setNameFilter(m_filter);
	if (d.exec()) {
		QStringList files = d.selectedFiles();
		m_last = d.directory();
		//_settings.setValue(_last_dir, m_last.path());
		std::set<std::string> tmp(m_files.begin(), m_files.end());
		for (int i = 0; i < files.size(); ++i)
			tmp.insert(files[i].toStdString());
		m_files.clear();
		m_files.assign(tmp.begin(), tmp.end());
	}
	updateFileList();
}

std::vector<std::string> FileList::files() {
	return m_files;
}

void FileList::setFiles(const std::vector<std::string> &files) {
	m_files = files;
	updateFileList();
}

FileList::~FileList() {
}

