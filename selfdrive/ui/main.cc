#include <csignal>
#include <execinfo.h>
#include <iostream>
#include <sys/resource.h>

#include <QApplication>
#include <QTranslator>

#include "system/hardware/hw.h"
#include "selfdrive/ui/qt/qt_window.h"
#include "selfdrive/ui/qt/util.h"
#include "selfdrive/ui/qt/window.h"

void crashHandler(int sig) {
  void *array[20];
  size_t size = backtrace(array, 20);
  std::cerr << "Error: signal " << sig << ":\n";
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

int main(int argc, char *argv[]) {
  signal(SIGABRT, crashHandler);
  signal(SIGSEGV, crashHandler);

  setpriority(PRIO_PROCESS, 0, -20);

  qInstallMessageHandler(swagLogMessageHandler);
  initApp(argc, argv);

  QTranslator translator;
  QString translation_file = QString::fromStdString(Params().get("LanguageSetting"));
  if (!translator.load(QString(":/%1").arg(translation_file)) && translation_file.length()) {
    qCritical() << "Failed to load translation file:" << translation_file;
  }

  QApplication a(argc, argv);
  a.installTranslator(&translator);

  MainWindow w;
  setMainWindow(&w);
  a.installEventFilter(&w);
  return a.exec();
}
