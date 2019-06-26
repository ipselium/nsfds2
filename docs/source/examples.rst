Examples
========

Here are some configuration examples. Some of these examples work with a :download:`geo file <examples/geo.py>`:

- :download:`Reference case <examples/reference.conf>`
- :download:`Harmonic source <examples/harmonic.conf>`
- :download:`Circle <examples/circle.conf>`
- :download:`Ground <examples/ground.conf>`
- :download:`Obstacles with moving boundary <examples/moving.conf>`
- :download:`Wave source <examples/wavefile.conf>` with :download:`a wavfile <examples/spam3.wav>`

To run the solver with one of these configuration::

   nsfds2 solve -c reference.conf
