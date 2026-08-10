/*
 * voidy-menu.js
 * Navbar dropdown behaviour for the multi-level menu (Bootstrap 5, no jQuery).
 *
 * Desktop (>= 992px):
 *   - submenus open on hover (pure CSS, see voidybootstrap-custom.css),
 *   - clicking a parent item navigates to its own page.
 * Mobile (collapsed navbar):
 *   - tapping a parent item unfolds its submenu inline instead of
 *     navigating (each submenu contains the parent page itself as its
 *     first entry, so every page stays reachable).
 */
(function () {
  'use strict';

  var DESKTOP = '(min-width: 992px)';

  function isDesktop() {
    return window.matchMedia(DESKTOP).matches;
  }

  function closeOthers(current) {
    document.querySelectorAll('#main-navbar .dropdown-menu.show').forEach(function (menu) {
      if (!menu.contains(current) && menu !== current) {
        menu.classList.remove('show');
        var toggle = menu.parentElement.querySelector(':scope > .dropdown-toggle');
        if (toggle) { toggle.setAttribute('aria-expanded', 'false'); }
      }
    });
  }

  document.querySelectorAll('#main-navbar .dropdown-toggle').forEach(function (toggle) {
    toggle.addEventListener('click', function (event) {
      if (isDesktop()) {
        // Desktop: hover already shows the submenu; a click follows the link.
        return;
      }
      // Mobile: toggle the submenu instead of navigating.
      event.preventDefault();
      event.stopPropagation();
      var parent = toggle.parentElement;
      var menu = parent.querySelector(':scope > .dropdown-menu');
      if (!menu) { return; }
      var shown = menu.classList.toggle('show');
      toggle.setAttribute('aria-expanded', shown ? 'true' : 'false');
      if (shown) { closeOthers(menu); }
    });
  });

  // Collapse all open submenus when the navbar collapse itself is closed.
  var collapse = document.getElementById('main-navbar-collapse');
  if (collapse) {
    collapse.addEventListener('hidden.bs.collapse', function () {
      collapse.querySelectorAll('.dropdown-menu.show').forEach(function (menu) {
        menu.classList.remove('show');
      });
    });
  }
})();
