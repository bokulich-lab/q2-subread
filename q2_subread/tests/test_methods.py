import unittest

from q2_subread.plugin_setup import plugin


class PluginSetupTests(unittest.TestCase):
    def test_plugin_metadata(self):
        self.assertEqual(plugin.name, 'subread')
        self.assertEqual(plugin.package, 'q2_subread')
        self.assertEqual(
            plugin.description,
            'QIIME 2 plugin for wrapping the Subread software package.',
        )

    def test_no_actions_registered_yet(self):
        self.assertEqual(len(plugin.methods), 0)
