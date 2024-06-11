import os
import sys
import subprocess
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QLineEdit, QPushButton, QFileDialog, QComboBox, QMessageBox, QVBoxLayout, QWidget
from PyQt5.QtGui import QPixmap, QColor, QPalette
from PyQt5.QtCore import Qt

class ViperGUI(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("VIPER Phylogeny")

        # Acrylic layout
        self.setWindowOpacity(0.95)  # Adjust alpha value for transparency

        # Expand the user's home directory in the logo path
        self.logo_path = os.path.expanduser("~/viperGUI/VIPER_ico.png")

        self.central_widget = QWidget(self)
        self.setCentralWidget(self.central_widget)

        self.folder_path = QLineEdit(self)
        self.num_threads = QComboBox(self)
        self.num_threads.setEditable(True)  # Make the combo box editable
        self.num_threads.setEditText("1")  # Default number of threads is set to 1

        # Virus Selection
        self.virus_choice = QComboBox(self)

        # UFBoostrap Selection
        self.ufbootstrap_choice = QLineEdit(self)
        self.ufbootstrap_choice.setPlaceholderText("1000")

        # Evolutionary Model Selection
        self.evolutionary_model_choice = QLineEdit(self)
        self.evolutionary_model_choice.setPlaceholderText("TEST")

        # Input Fasta Sequences
        self.input_fasta_sequences = QLineEdit(self)

        # Metadata TSV File
        self.metadata_tsv_file = QLineEdit(self)

        # Analysis Name
        self.analysis_name = QLineEdit(self)

        # Create and set the logo
        self.logo = QLabel(self)
        self.logo.setAlignment(Qt.AlignTop)

        # Create and set up widgets
        self.create_widgets()

    def create_widgets(self):
        layout = QVBoxLayout(self.centralWidget())

        # Logo
        self.update_logo_size()  # Update logo size initially
        layout.addWidget(self.logo)

        # Folder Selection
        label_folder = QLabel("Select Folder:", self)
        layout.addWidget(label_folder)
        layout.addWidget(self.folder_path)
        button_browse = QPushButton("Browse", self)
        button_browse.clicked.connect(self.browse_folder)
        layout.addWidget(button_browse)

        # Virus Selection
        label_virus = QLabel("Select Virus:", self)
        layout.addWidget(label_virus)
        layout.addWidget(self.virus_choice)
        viruses = [
            "SARS-CoV-2", "DENV-1", "DENV-2", "DENV-3", "DENV-4",
            "Influenza A/H1N1 - HA", "Influenza A/H1N1 - NA",
            "Influenza A/H3N2 - HA", "Influenza A/H3N2 - NA",
            "Influenza B/Victoria - HA", "Influenza B/Victoria - NA",
            "Influenza B/Yamagata- HA", "Influenza B/Yamagata - NA"
        ]
        self.virus_choice.addItems(viruses)

        # UFBoostrap Selection
        label_ufbootstrap = QLabel("Select UFBoostrap:", self)
        layout.addWidget(label_ufbootstrap)
        layout.addWidget(self.ufbootstrap_choice)

        # Evolutionary Model Selection
        label_evolutionary_model = QLabel("Select Evolutionary Model:", self)
        layout.addWidget(label_evolutionary_model)
        layout.addWidget(self.evolutionary_model_choice)

        # Input Fasta Sequences
        label_sequences = QLabel("Input Fasta Sequences:", self)
        layout.addWidget(label_sequences)
        layout.addWidget(self.input_fasta_sequences)
        button_browse_sequences = QPushButton("Browse", self)
        button_browse_sequences.clicked.connect(self.browse_sequences)
        layout.addWidget(button_browse_sequences)

        # Metadata TSV File
        label_metadata_tsv = QLabel("Metadata TSV File (Optional):", self)
        layout.addWidget(label_metadata_tsv)
        layout.addWidget(self.metadata_tsv_file)
        button_browse_metadata = QPushButton("Browse", self)
        button_browse_metadata.clicked.connect(self.browse_metadata)
        layout.addWidget(button_browse_metadata)

        # Analysis Name
        label_analysis_name = QLabel("Analysis Name:", self)
        layout.addWidget(label_analysis_name)
        layout.addWidget(self.analysis_name)

        # Number of Threads Selection
        label_threads = QLabel("Number of Threads to Use:", self)
        layout.addWidget(label_threads)
        layout.addWidget(self.num_threads)
        available_threads = os.cpu_count() or 1
        thread_values = [str(i) for i in range(1, available_threads + 1)]
        self.num_threads.addItems(thread_values)
        self.num_threads.currentIndexChanged.connect(self.on_threads_changed)

        # Run Button
        button_run = QPushButton("Run Phylogenetic analysis", self)
        button_run.clicked.connect(self.run_pipeline)
        layout.addWidget(button_run)

        # Update Mutation Tables Button
        button_update_mutations = QPushButton("Update Mutation Tables", self)
        button_update_mutations.clicked.connect(self.update_mutation_tables)
        layout.addWidget(button_update_mutations)

        # Apply style to adjust font size
        self.setStyleSheet(
            """
            QLabel {
                font-size: 14px;
            }
            QPushButton {
                font-size: 14px;
                padding: 8px;
            }
            """
        )

        # Apply Fusion style
        app.setStyle("Fusion")
        dark_palette = QPalette()
        dark_palette.setColor(QPalette.Window, QColor(41, 41, 41))
        dark_palette.setColor(QPalette.WindowText, Qt.white)
        dark_palette.setColor(QPalette.Base, QColor(25, 25, 25))
        dark_palette.setColor(QPalette.AlternateBase, QColor(41, 41, 41))
        dark_palette.setColor(QPalette.ToolTipBase, Qt.white)
        dark_palette.setColor(QPalette.ToolTipText, Qt.white)
        dark_palette.setColor(QPalette.Text, Qt.white)
        dark_palette.setColor(QPalette.Button, QColor(41, 41, 41))
        dark_palette.setColor(QPalette.ButtonText, Qt.white)
        dark_palette.setColor(QPalette.BrightText, Qt.red)
        dark_palette.setColor(QPalette.Link, QColor(42, 130, 218))
        dark_palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
        dark_palette.setColor(QPalette.HighlightedText, Qt.black)
        app.setPalette(dark_palette)

    def browse_folder(self):
        folder_selected = QFileDialog.getExistingDirectory(self, "Select Folder")
        self.folder_path.setText(folder_selected)

    def browse_sequences(self):
        file_selected, _ = QFileDialog.getOpenFileName(self, "Select Fasta Sequences File", filter="Fasta Files (*.fasta *.fa)")
        self.input_fasta_sequences.setText(file_selected)

    def browse_metadata(self):
        file_selected, _ = QFileDialog.getOpenFileName(self, "Select Metadata TSV File", filter="TSV Files (*.tsv)")
        self.metadata_tsv_file.setText(file_selected)

    def run_pipeline(self):
        folder = self.folder_path.text()
        num_threads = self.num_threads.currentText()
        virus = self.virus_choice.currentText()
        ufbootstrap = self.ufbootstrap_choice.text() or "1000"
        evolutionary_model = self.evolutionary_model_choice.text() or "TEST"
        input_fasta_sequences = self.input_fasta_sequences.text()
        metadata_tsv_file = self.metadata_tsv_file.text() or "null"
        analysis_name = self.analysis_name.text()

        if not folder or not input_fasta_sequences or not analysis_name:
            QMessageBox.critical(self, "Error", "Please fill in all required fields.")
            return

        try:
            # Change the working directory to the selected folder
            os.chdir(folder)
            # Run the pipeline script with the specified parameters
            command = f"execPhylo.sh -i {input_fasta_sequences} -v {self.get_virus_alias(virus)} -j {analysis_name} -n {num_threads} -metadata {metadata_tsv_file} -b {ufbootstrap} -m {evolutionary_model}"
            print(f"Running command: {command}")
            subprocess.run(command, check=True, shell=True)

            QMessageBox.information(self, "Executed", "VIPER executed!")
        except subprocess.CalledProcessError:
            QMessageBox.critical(self, "Error", "An error occurred while trying to run VIPER. Please, check the screen log.")

    def update_mutation_tables(self):
        pipeline_path = os.environ.get("PIPELINE")
        if not pipeline_path:
            QMessageBox.critical(self, "Error", "Pipeline path environment variable (PIPELINE) not set.")
            return

        script_path = os.path.join(pipeline_path, "update_database", "updateMutationTables.sh")

        try:
            # Run the updateMutationTables.sh script
            command = f"sh {script_path}"
            print(f"Running command: {command}")
            subprocess.run(command, check=True, shell=True)
            QMessageBox.information(self, "Success", "Mutation tables updated successfully!")
        except subprocess.CalledProcessError:
            QMessageBox.critical(self, "Error", "An error occurred while trying to update mutation tables. Please, check the screen log.")

    def on_threads_changed(self, index):
        selected_threads = self.num_threads.itemText(index)
        print(f"Selected Threads: {selected_threads}")

    def get_virus_alias(self, virus):
        virus_aliases = {
            "Influenza A/H1N1 - HA": "FLUA_H1",
            "Influenza A/H1N1 - NA": "FLUA_N1",
            "Influenza A/H3N2 - HA": "FLUB_H3",
            "Influenza A/H3N2 - NA": "FLUB_N2",
            "Influenza B/Victoria - HA": "FLUB_VIC_HA",
            "Influenza B/Victoria - NA": "FLUB_VIC_NA",
            "Influenza B/Yamagata- HA": "FLUB_YAM_HA",
            "Influenza B/Yamagata - NA": "FLUB_YAM_NA"
        }
        return virus_aliases.get(virus, virus)

    def resizeEvent(self, event):
        # Called when the window is resized
        self.update_logo_size()

    def update_logo_size(self):
        # Update logo size based on window width
        logo_width = min(self.width() - 20, 500)  # Limit maximum width to 500 pixels
        self.logo.setPixmap(QPixmap(self.logo_path).scaledToWidth(logo_width, Qt.SmoothTransformation))

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = ViperGUI()
    window.setGeometry(100, 100, 520, 500)  # Adjusted window height to accommodate new options
    window.show()
    sys.exit(app.exec_())
