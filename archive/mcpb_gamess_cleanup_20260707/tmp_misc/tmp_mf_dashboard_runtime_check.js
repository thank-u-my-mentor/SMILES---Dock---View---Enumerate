const fs = require("fs");
const vm = require("vm");

const scriptPath = "/mnt/e/Codex/tmp_mf_dashboard_script.js";
const source = fs.readFileSync(scriptPath, "utf8");
const elements = {};
function element(id) {
  if (!elements[id]) {
    elements[id] = {
      id,
      value: "",
      innerHTML: "",
      textContent: "",
      addEventListener() {},
    };
  }
  return elements[id];
}

const document = {
  getElementById: element,
  querySelectorAll() {
    return [];
  },
};
const localStorage = {
  getItem() {
    return "";
  },
  setItem() {},
};
const context = {
  document,
  localStorage,
  setTimeout(fn) {
    fn();
  },
  console,
};

vm.createContext(context);
vm.runInContext(
  source +
    "\nactiveTab='identity'; renderDetail();" +
    "\nif(!document.getElementById('list').innerHTML.includes('SpPPase')) throw new Error('gene label missing from list');" +
    "\nif(!document.getElementById('detail').innerHTML.includes('identity-heatmap')) throw new Error('identity heatmap missing');" +
    "\nif(!document.getElementById('detail').innerHTML.includes('RnARG')) throw new Error('matrix gene label missing');",
  context,
);
console.log("runtime check passed");
