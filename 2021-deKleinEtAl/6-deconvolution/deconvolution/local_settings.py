"""
File:         local_settings.py
Created:      2020/10/08
Last Changed:
Author:       M.Vochteloo

Copyright (C) 2020 M.Vochteloo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License can be found in the LICENSE file in the
root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
"""

# Standard imports.
import json
import os

# Third party imports.

# Local application imports.


class LocalSettings:
    class __LocalSettings:
        def __init__(self, project_root, settings_file):
            self.settings = None
            self.settings_path = os.path.join(project_root,
                                              'settings',
                                              settings_file + ".json")
            self.load_settings()

        def load_settings(self):
            settings_content = open(self.settings_path)
            settings = json.loads(settings_content.read())
            settings_content.close()
            self.settings = settings

        def is_settings_loaded(self):
            return self.settings is not None

        def get_setting(self, setting_key):
            if self.settings is not None and setting_key in self.settings:
                value = self.settings[setting_key]
            else:
                value = None

            return value

        def get_all_settings(self):
            return self.settings

        def __str__(self):
            return repr(self)

    instance = None

    def __new__(cls, project_root, settings_file):
        if not LocalSettings.instance:
            LocalSettings.instance = LocalSettings.__LocalSettings(project_root,
                                                                   settings_file)
        return LocalSettings.instance

    def __init__(self, project_root, settings_file):
        if not LocalSettings.instance:
            LocalSettings.instance = LocalSettings.__LocalSettings(project_root,
                                                                   settings_file)

    def __getattr__(self, name):
        return getattr(self.instance, name)
