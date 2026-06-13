import sys
import duckdb
import pyperclip
from textual import on, events
from textual.app import App, ComposeResult
from textual.containers import Horizontal, Vertical, HorizontalScroll
from textual.widgets import Input, DataTable, Label, Button
from textual.reactive import reactive

MAX_CELL_LENGTH = 20

class VimInput(Input):
    def on_focus(self, event: events.Focus) -> None:
        self.app.last_focused_input = self
        try:
            self.scroll_visible(center=True, animate=True)
        except TypeError:
            self.scroll_visible(animate=True)

    def on_key(self, event: events.Key) -> None:
        if self.app.current_mode in ("normal", "command"):
            event.prevent_default()
            
            # Only process navigation if we are strictly in normal mode
            if self.app.current_mode == "normal":
                if event.key in ("h", "j", "k", "l", "i", "c"):
                    event.stop() 
                    
                    if event.key == "i":
                        self.app.current_mode = "insert"
                    elif event.key == "j":
                        btn_prev = self.app.query_one("#btn_prev")
                        btn_next = self.app.query_one("#btn_next")
                        if not btn_prev.disabled:
                            btn_prev.focus()
                        elif not btn_next.disabled:
                            btn_next.focus()
                        else:
                            self.app.query_one(VimDataTable).focus()
                    elif event.key in ("h", "l"):
                        inputs = list(self.app.query(VimInput))
                        idx = inputs.index(self)
                        if event.key == "h":
                            target_idx = (idx - 1) % len(inputs)
                            inputs[target_idx].focus()
                        elif event.key == "l":
                            target_idx = (idx + 1) % len(inputs)
                            inputs[target_idx].focus()

class VimButton(Button):
    def on_key(self, event: events.Key) -> None:
        if self.app.current_mode == "normal":
            if event.key in ("h", "j", "k", "l"):
                event.prevent_default()
                event.stop()  # Prevents the key from cascading
                
                if event.key == "j":
                    self.app.query_one(VimDataTable).focus()
                elif event.key == "k":
                    last_input = getattr(self.app, "last_focused_input", None)
                    if last_input:
                        last_input.focus()
                    else:
                        self.app.query(VimInput).first().focus()
                elif event.key in ("h", "l"):
                    btn_prev = self.app.query_one("#btn_prev")
                    btn_next = self.app.query_one("#btn_next")
                    
                    if self.id == "btn_prev" and event.key == "l" and not btn_next.disabled:
                        btn_next.focus()
                    elif self.id == "btn_next" and event.key == "h" and not btn_prev.disabled:
                        btn_prev.focus()

class VimDataTable(DataTable):
    BINDINGS = [
        ("h", "cursor_left", "Left"),
        ("j", "cursor_down", "Down"),
        ("k", "custom_up", "Up"),
        ("l", "cursor_right", "Right"),
        ("i", "custom_insert", "Insert Mode"),
    ]

    def on_focus(self, event: events.Focus) -> None:
        self.show_cursor = True

    def on_blur(self, event: events.Blur) -> None:
        self.show_cursor = False

    def action_custom_up(self) -> None:
        if self.cursor_coordinate.row == 0:
            btn_prev = self.app.query_one("#btn_prev")
            btn_next = self.app.query_one("#btn_next")
            
            if not btn_prev.disabled:
                btn_prev.focus()
            elif not btn_next.disabled:
                btn_next.focus()
            else:
                last_input = getattr(self.app, "last_focused_input", None)
                if last_input:
                    last_input.focus()
                else:
                    self.app.query(VimInput).first().focus()
        else:
            self.action_cursor_up()

    def action_custom_insert(self) -> None:
        self.app.current_mode = "insert"
        last_input = getattr(self.app, "last_focused_input", None)
        if last_input:
            last_input.focus()
        else:
            self.app.query(VimInput).first().focus()

class ColumnSearchApp(App):
    CSS = """
        /* Tokyo Night Theme Variables - Must be at root */
        $background: #1a1b26;
        $surface: #24283b;
        $surface-light: #292e42;
        $panel: #24283b;

        $primary: #7aa2f7;
        $primary-light: #7dcfff;
        $primary-dark: #3d59a1;

        $secondary: #bb9af7;
        $accent: #ff9e64;

        $success: #9ece6a;
        $warning: #e0af68;
        $error: #f7768e;

        $text: #c0caf5;
        $text-muted: #565f89;

        Screen {
          background: $background;
        }

        #search-container {
          layout: horizontal;
          align: left middle;
          height: auto;
          padding: 1;
          border-bottom: solid $primary-dark;
          background: $surface;
        }

        .filter-box {
          width: 30;
          height: auto;
          margin: 0 1 2 1;
        }

        .col-label {
          margin-left: 1;
          text-style: bold;
          color: $primary-light;
        }

        Input {
          background: $surface-light;
          border: tall $surface-light;
          color: $text;
        }

        Input:focus {
          border: tall $primary;
        }

        #pagination-container {
          layout: horizontal;
          align: center middle;
          height: auto;
          padding: 1;
        }

        #pagination-container VimButton {
          min-width: 10;
          margin: 0 2;
        }

        #result-count {
          content-align: center middle;
          padding: 1 2;
          text-style: bold;
          color: $success;
        }

        #mode-indicator {
          content-align: center middle;
          text-style: bold;
          color: $background;
          padding: 0 2;
        }

        #mode-indicator.mode-normal { background: $primary; }
        #mode-indicator.mode-insert { background: $success; }
        #mode-indicator.mode-command { background: $warning; }

        VimDataTable {
          width: 1fr;
          height: 1fr;
          background: $background;
          color: $text;
        }

        VimDataTable > .datatable--header {
          background: $surface;
          color: $primary;
          text-style: bold;
        }

        VimDataTable > .datatable--cursor {
          background: $primary-dark;
          color: $text;
        }

        /* Custom Status Bar CSS */
        #status-bar {
          dock: bottom;
          height: 1;
          background: $surface;
          color: $text;
          layout: horizontal;
        }

        #mode-indicator {
          content-align: center middle;
          text-style: bold;
          color: $background;
          background: $success;
          padding: 0 2;
        }

        #status-hint {
          content-align: center middle;
          padding: 0 2;
          color: $text-muted;
        }
    """ 
    BINDINGS = [
        ("escape", "switch_normal", "Normal Mode"),
        ("c", "copy_cell", "Copy Original Cell Value")
    ]

    current_mode = reactive("normal")
    command_sequence = ""

    def watch_current_mode(self, old_mode: str, new_mode: str) -> None:
        try:
            indicator = self.query_one("#mode-indicator", Label)
            indicator.update(f" {new_mode.upper()} ")
            
            if new_mode == "normal":
                indicator.styles.background = "#7aa2f7"
                for inp in self.query(VimInput):
                    inp.cursor_blink = False
            elif new_mode == "command":
                indicator.styles.background = "#e0af68" # Yellow/Warning color
                for inp in self.query(VimInput):
                    inp.cursor_blink = False
            else:
                indicator.styles.background = "#9ece6a"
                for inp in self.query(VimInput):
                    inp.cursor_blink = True
                    if inp.has_focus:
                        inp.refresh()
        except Exception:
            pass

    def on_key(self, event: events.Key) -> None:
        # Handle state transitions that bubble up
        if self.current_mode == "normal":
            if event.character == ":":
                self.current_mode = "command"
        elif self.current_mode == "command":
            if event.character == "q":
                self.exit()
            else:
                self.current_mode = "normal"

    def __init__(self, db_path: str):
        super().__init__()
        self.db_path = db_path
        self.con = duckdb.connect(database=self.db_path, read_only=True)
        
        tables = self.con.execute("SHOW TABLES").fetchall()
        self.table_name = tables[0][0] if tables else "data"
        
        self.columns = [col[0] for col in self.con.execute(f"DESCRIBE {self.table_name}").fetchall()]
        self.current_data = []
        
        self.limit = 25
        self.offset = 0
        self.total_count = 0

    def compose(self) -> ComposeResult:
        with Vertical():
            with HorizontalScroll(id="search-container"):
                for col in self.columns:
                    with Vertical(classes="filter-box"):
                        yield Label(col, classes="col-label")
                        yield VimInput(id=f"input_{col}")
            
            with Horizontal(id="pagination-container"):
                yield VimButton("Prev", id="btn_prev", disabled=True)
                yield Label("Results: 0", id="result-count")
                yield VimButton("Next", id="btn_next", disabled=True)
            
            yield VimDataTable(id="table")
            
        with Horizontal(id="status-bar"):
            yield Label(" NORMAL ", id="mode-indicator", classes="mode-normal")
            yield Label(" ESC: Normal | i: Insert | :q: Quit | hjkl: Navigate", id="status-hint")

    def on_mount(self) -> None:
        table = self.query_one(VimDataTable)
        table.add_columns(*self.columns)
        table.cursor_type = "cell" 
        table.show_cursor = False  # Start with the cursor hidden
        self.update_table()
       
        for inp in self.query(VimInput):
            inp.cursor_blink = False

        self.query(VimInput).first().focus()

    def watch_current_mode(self, old_mode: str, new_mode: str) -> None:
        try:
            indicator = self.query_one("#mode-indicator", Label)
            indicator.update(f" {new_mode.upper()} ")
            
            # Dynamically swap the CSS classes
            indicator.remove_class(f"mode-{old_mode}")
            indicator.add_class(f"mode-{new_mode}")
            
            # Handle cursor blink
            if new_mode in ("normal", "command"):
                for inp in self.query(VimInput):
                    inp.cursor_blink = False
            else:
                for inp in self.query(VimInput):
                    inp.cursor_blink = True
                    if inp.has_focus:
                        inp.refresh()
        except Exception:
            pass
 
    def action_switch_normal(self) -> None:
        self.current_mode = "normal"

    @on(Input.Changed)
    def handle_input_change(self) -> None:
        self.offset = 0 
        self.update_table()

    @on(Button.Pressed, "#btn_prev")
    def handle_prev_page(self) -> None:
        self.offset = max(0, self.offset - self.limit)
        self.update_table()

    @on(Button.Pressed, "#btn_next")
    def handle_next_page(self) -> None:
        if self.offset + self.limit < self.total_count:
            self.offset += self.limit
            self.update_table()

    def action_copy_cell(self) -> None:
        table = self.query_one(VimDataTable)
        if not self.current_data:
            return
            
        try:
            row_idx = table.cursor_coordinate.row
            col_idx = table.cursor_coordinate.column
            
            original_value = self.current_data[row_idx][col_idx]
            
            pyperclip.copy(str(original_value))
            self.notify("Copied to clipboard!") 
        except Exception:
            pass

    def update_table(self) -> None:
        table = self.query_one(VimDataTable)
        table.clear()

        conditions = []
        for col in self.columns:
            val = self.query_one(f"#input_{col}", VimInput).value.strip()
            if val:
                safe_val = val.replace("'", "''")
                conditions.append(f"{col}::VARCHAR ILIKE '%{safe_val}%'")

        where_clause = ""
        if conditions:
            where_clause = " WHERE " + " AND ".join(conditions)

        try:
            count_query = f"SELECT COUNT(*) FROM {self.table_name}{where_clause}"
            self.total_count = self.con.execute(count_query).fetchone()[0]

            data_query = f"SELECT * FROM {self.table_name}{where_clause} LIMIT {self.limit} OFFSET {self.offset}"
            self.current_data = self.con.execute(data_query).fetchall()
            
            current_page = (self.offset // self.limit) + 1
            total_pages = max(1, (self.total_count + self.limit - 1) // self.limit)
            
            self.query_one("#result-count", Label).update(
                f"Page {current_page} of {total_pages} | Total: {self.total_count}"
            )
            
            self.query_one("#btn_prev", VimButton).disabled = self.offset == 0
            self.query_one("#btn_next", VimButton).disabled = current_page >= total_pages

            display_data = []
            for row in self.current_data:
                display_row = [
                    f"{str(x)[:MAX_CELL_LENGTH]}..." if len(str(x)) > MAX_CELL_LENGTH else str(x)
                    for x in row
                ]
                display_data.append(display_row)

            table.add_rows(display_data)
        except Exception:
            self.query_one("#result-count", Label).update("Error querying database")

def main():
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else "data/hypersurfaces/hypersurfaces.db"
    app = ColumnSearchApp(path)
    app.run()

if __name__ == "__main__":
    main()
