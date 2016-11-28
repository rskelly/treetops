#include "fileselectorpanel.h"
#include "fileselectorpanelplugin.h"

#include <QtPlugin>

FileSelectorPanelPlugin::FileSelectorPanelPlugin(QObject *parent)
    : QObject(parent)
{
    m_initialized = false;
}

void FileSelectorPanelPlugin::initialize(QDesignerFormEditorInterface * /* core */)
{
    if (m_initialized)
        return;

    // Add extension registrations, etc. here

    m_initialized = true;
}

bool FileSelectorPanelPlugin::isInitialized() const
{
    return m_initialized;
}

QWidget *FileSelectorPanelPlugin::createWidget(QWidget *parent)
{
    return new FileSelectorPanel(parent);
}

QString FileSelectorPanelPlugin::name() const
{
    return QLatin1String("FileSelectorPanel");
}

QString FileSelectorPanelPlugin::group() const
{
    return QLatin1String("geotools");
}

QIcon FileSelectorPanelPlugin::icon() const
{
    return QIcon();
}

QString FileSelectorPanelPlugin::toolTip() const
{
    return QLatin1String("");
}

QString FileSelectorPanelPlugin::whatsThis() const
{
    return QLatin1String("");
}

bool FileSelectorPanelPlugin::isContainer() const
{
    return false;
}

QString FileSelectorPanelPlugin::domXml() const
{
    return QLatin1String("<widget class=\"FileSelectorPanel\" name=\"fileSelectorPanel\">\n</widget>\n");
}

QString FileSelectorPanelPlugin::includeFile() const
{
    return QLatin1String("fileselectorpanel.h");
}
#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(fileselectorpanelplugin, FileSelectorPanelPlugin)
#endif // QT_VERSION < 0x050000
