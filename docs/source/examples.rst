Examples
========

Here are some configuration examples. Some of these examples work with a :download:`geo file <examples/geo.py>`:

- :download:`Reference case <examples/reference.conf>`
- :download:`Letter A <examples/custom.conf>`
  (:download:`geo file <examples/geo.py>` needed)
- :download:`Harmonic source <examples/harmonic.conf>`
- :download:`Circle <examples/circle.conf>`
- :download:`Ground <examples/ground.conf>`
- :download:`Periodic <examples/periodic.conf>`
- :download:`Section change <examples/section.conf>`
  (:download:`geo file <examples/geo.py>` needed)
- :download:`Obstacles with moving boundary <examples/moving.conf>`
  (:download:`geo file <examples/geo.py>` needed)
- :download:`Wave source <examples/wavefile.conf>` with
  :download:`Wavfile example <examples/spam3.wav>`
  (:download:`geo file <examples/geo.py>` needed)


To run the solver with one of these configuration::

   nsfds2 solve -c reference.conf
