use std::fs;
use std::sync::Mutex;

use tauri::{AppHandle, Emitter, Manager, State};
use tauri_plugin_dialog::DialogExt;
use tauri_plugin_shell::process::CommandEvent;
use tauri_plugin_shell::ShellExt;

struct ProcessState {
    child: Mutex<Option<tauri_plugin_shell::process::CommandChild>>,
}

#[derive(Clone, serde::Serialize)]
struct StatusPayload {
    phase: String,
    progress: i32,
    message: String,
}

#[derive(Clone, serde::Serialize)]
struct LogPayload {
    line: String,
    stream: String,
}

#[tauri::command]
fn read_settings_file(path: String) -> Result<String, String> {
    fs::read_to_string(&path).map_err(|e| format!("Failed to read {path}: {e}"))
}

#[tauri::command]
fn write_settings_file(path: String, content: String) -> Result<(), String> {
    if let Some(parent) = std::path::Path::new(&path).parent() {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    fs::write(&path, content).map_err(|e| format!("Failed to write {path}: {e}"))
}

#[tauri::command]
async fn pick_open_file(app: AppHandle, title: String) -> Result<Option<String>, String> {
    let path = app
        .dialog()
        .file()
        .set_title(&title)
        .blocking_pick_file();
    Ok(path.map(|p| p.to_string()))
}

#[tauri::command]
async fn pick_save_file(app: AppHandle, title: String) -> Result<Option<String>, String> {
    let path = app
        .dialog()
        .file()
        .set_title(&title)
        .blocking_save_file();
    Ok(path.map(|p| p.to_string()))
}

#[tauri::command]
async fn cancel_processing(state: State<'_, ProcessState>) -> Result<(), String> {
    let mut guard = state.child.lock().map_err(|e| e.to_string())?;
    if let Some(child) = guard.take() {
        child.kill().map_err(|e| e.to_string())?;
    }
    Ok(())
}

#[tauri::command]
async fn run_processing(
    app: AppHandle,
    config_path: String,
    state: State<'_, ProcessState>,
) -> Result<(), String> {
    {
        let guard = state.child.lock().map_err(|e| e.to_string())?;
        if guard.is_some() {
            return Err("Processing is already running.".into());
        }
    }

    let sidecar = app
        .shell()
        .sidecar("treetops-cli")
        .map_err(|e| format!("Sidecar not found. Run npm run prepare-sidecar after building treetops-cli: {e}"))?;

    let (mut rx, child) = sidecar
        .args(["-c", &config_path])
        .spawn()
        .map_err(|e| e.to_string())?;

    {
        let mut guard = state.child.lock().map_err(|e| e.to_string())?;
        *guard = Some(child);
    }

    let app_handle = app.clone();
    tauri::async_runtime::spawn(async move {
        let mut exit_code: Option<i32> = None;

        while let Some(event) = rx.recv().await {
            match event {
                CommandEvent::Stdout(bytes) => {
                    let line = String::from_utf8_lossy(&bytes).trim_end().to_string();
                    if let Some(json) = line.strip_prefix("@status:") {
                        if let Ok(value) = serde_json::from_str::<serde_json::Value>(json) {
                            let payload = StatusPayload {
                                phase: value
                                    .get("phase")
                                    .and_then(|v| v.as_str())
                                    .unwrap_or("unknown")
                                    .to_string(),
                                progress: value
                                    .get("progress")
                                    .and_then(|v| v.as_i64())
                                    .unwrap_or(0) as i32,
                                message: value
                                    .get("message")
                                    .and_then(|v| v.as_str())
                                    .unwrap_or("")
                                    .to_string(),
                            };
                            let _ = app_handle.emit("status", payload);
                        }
                    } else if !line.is_empty() {
                        let _ = app_handle.emit(
                            "log",
                            LogPayload {
                                line,
                                stream: "stdout".into(),
                            },
                        );
                    }
                }
                CommandEvent::Stderr(bytes) => {
                    let line = String::from_utf8_lossy(&bytes).trim_end().to_string();
                    if !line.is_empty() {
                        let _ = app_handle.emit(
                            "log",
                            LogPayload {
                                line,
                                stream: "stderr".into(),
                            },
                        );
                    }
                }
                CommandEvent::Terminated(payload) => {
                    exit_code = payload.code;
                }
                _ => {}
            }
        }

        if let Some(state) = app_handle.try_state::<ProcessState>() {
            if let Ok(mut guard) = state.child.lock() {
                *guard = None;
            }
        }

        let success = exit_code == Some(0);
        let message = if success {
            "Processing complete.".into()
        } else {
            format!("Processing failed (exit code {:?}).", exit_code)
        };

        let _ = app_handle.emit(
            "status",
            StatusPayload {
                phase: if success { "done" } else { "error" }.into(),
                progress: if success { 100 } else { 0 },
                message,
            },
        );
        let _ = app_handle.emit("processing-finished", success);
    });

    Ok(())
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_shell::init())
        .manage(ProcessState {
            child: Mutex::new(None),
        })
        .invoke_handler(tauri::generate_handler![
            read_settings_file,
            write_settings_file,
            pick_open_file,
            pick_save_file,
            run_processing,
            cancel_processing
        ])
        .run(tauri::generate_context!())
        .expect("error while running treetops application");
}
