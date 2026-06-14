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
    """Minimal subclass to prevent text entry when not in insert mode."""
    
    def on_key(self, event: events.Key) -> None:
        if self.app.current_mode != "insert":
            # Block standard printable characters
            if event.is_printable:
                event.prevent_default()

    def action_delete_left(self) -> None:
        """Fired by the backspace key."""
        if self.app.current_mode == "insert":
            super().action_delete_left()

    def action_delete_right(self) -> None:
        """Fired by the delete key."""
        if self.app.current_mode == "insert":
            super().action_delete_right()

class ColumnSearchApp(App):
    CSS = """
        /* Tokyo Night Theme Variables */
        $background: #1a1b26;
        $surface: #24283b;
        $surface-light: #292e42;
        $panel: #24283b;
        $primary: #7aa2f7;
        $primary-light: #7dcfff;
        $primary-dark: #3d59a1;
        $success: #9ece6a;
        $warning: #e0af68;
        $text: #c0caf5;
        $text-muted: #565f89;

        Screen { background: $background; }

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
        Input:focus { border: tall $primary; }

        #pagination-container {
          layout: horizontal;
          align: center middle;
          height: auto;
          padding: 1;
        }

        #pagination-container Button {
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

        DataTable {
          width: 1fr;
          height: 1fr;
          background: $background;
          color: $text;
        }

        DataTable > .datatable--header {
          background: $surface;
          color: $primary;
          text-style: bold;
        }

        /* Hide the cursor when the table is unfocused */
        DataTable > .datatable--cursor {
          background: transparent;
          color: $text;
        }

        /* Show the cursor only when the table has focus */
        DataTable:focus > .datatable--cursor {
          background: $primary-dark;
          color: $text;
        }

        #status-bar {
          dock: bottom;
          height: 1;
          background: $surface;
          color: $text;
          layout: horizontal;
        }

        #status-hint {
          content-align: center middle;
          padding: 0 2;
          color: $text-muted;
        }
    """ 
    BINDINGS = [
        ("escape", "switch_normal", "Normal Mode"),
    ]

    current_mode = reactive("normal")

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
                yield Button("Prev", id="btn_prev", disabled=True)
                yield Label("Results: 0", id="result-count")
                yield Button("Next", id="btn_next", disabled=True)
            
            yield DataTable(id="table")
            
        with Horizontal(id="status-bar"):
            yield Label(" NORMAL ", id="mode-indicator", classes="mode-normal")
            yield Label(" ESC: Normal | i: Insert | :q: Quit | hjkl: Navigate | c: Copy", id="status-hint")

    def on_mount(self) -> None:
        table = self.query_one(DataTable)
        table.add_columns(*self.columns)
        table.cursor_type = "cell" 
        self.update_table()
        
        for inp in self.query(VimInput):
            inp.cursor_blink = False

        if inputs := self.query(VimInput):
            inputs.first().focus()

    def watch_current_mode(self, old_mode: str, new_mode: str) -> None:
        try:
            indicator = self.query_one("#mode-indicator", Label)
            indicator.update(f" {new_mode.upper()} ")
            indicator.remove_class(f"mode-{old_mode}")
            indicator.add_class(f"mode-{new_mode}")
            
            for inp in self.query(VimInput):
                inp.cursor_blink = (new_mode == "insert")
                if inp.has_focus:
                    inp.refresh()
        except Exception:
            pass

    def action_switch_normal(self) -> None:
        self.current_mode = "normal"

    def on_key(self, event: events.Key) -> None:
        if self.current_mode == "normal":
            if event.character == ":":
                self.current_mode = "command"
                return

            if event.key in ("h", "j", "k", "l", "i", "c"):
                if event.key == "i":
                    self.current_mode = "insert"
                    if not isinstance(self.focused, VimInput):
                        if inputs := self.query(VimInput):
                            inputs.first().focus()
                elif event.key == "c":
                    self.action_copy_cell()
                else:
                    self.handle_vim_navigation(event.key)
                    
        elif self.current_mode == "command":
            if event.character == "q":
                self.exit()
            elif event.key not in ("enter", "colon"):
                self.current_mode = "normal"

    def handle_vim_navigation(self, key: str) -> None:
        focused = self.focused
        
        if isinstance(focused, DataTable):
            try:
                if key == "j":
                    focused.action_cursor_down()
                elif key == "k":
                    if focused.cursor_coordinate.row <= 0:
                        btn_prev = self.query_one("#btn_prev", Button)
                        btn_next = self.query_one("#btn_next", Button)
                        if not btn_prev.disabled:
                            btn_prev.focus()
                        elif not btn_next.disabled:
                            btn_next.focus()
                        elif inputs := self.query(VimInput):
                            inputs.first().focus()
                    else:
                        focused.action_cursor_up()
                elif key == "h":
                    focused.action_cursor_left()
                elif key == "l":
                    focused.action_cursor_right()
            except Exception:
                pass

        elif isinstance(focused, Button):
            if key == "j":
                self.query_one(DataTable).focus()
            elif key == "k":
                if inputs := self.query(VimInput):
                    inputs.first().focus()
            elif key == "h" and focused.id == "btn_next":
                btn_prev = self.query_one("#btn_prev", Button)
                if not btn_prev.disabled:
                    btn_prev.focus()
            elif key == "l" and focused.id == "btn_prev":
                btn_next = self.query_one("#btn_next", Button)
                if not btn_next.disabled:
                    btn_next.focus()

        elif isinstance(focused, VimInput):
            if key == "j":
                btn_prev = self.query_one("#btn_prev", Button)
                btn_next = self.query_one("#btn_next", Button)
                if not btn_prev.disabled:
                    btn_prev.focus()
                elif not btn_next.disabled:
                    btn_next.focus()
                else:
                    self.query_one(DataTable).focus()
            elif key in ("h", "l"):
                inputs = list(self.query(VimInput))
                if inputs:
                    try:
                        idx = inputs.index(focused)
                    except ValueError:
                        idx = 0
                    
                    target_idx = (idx - 1) if key == "h" else (idx + 1)
                    inputs[target_idx % len(inputs)].focus()

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
        table = self.query_one(DataTable)
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
        table = self.query_one(DataTable)
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
            
            self.query_one("#btn_prev", Button).disabled = self.offset == 0
            self.query_one("#btn_next", Button).disabled = current_page >= total_pages

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
    path = sys.argv[1] if len(sys.argv) > 1 else "data/hypersurfaces/hypersurfaces.db"
    app = ColumnSearchApp(path)
    app.run()

if __name__ == "__main__":
    main()
