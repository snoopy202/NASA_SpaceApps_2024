import * as THREE from 'three';

// Initialize scene, camera, and renderer
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
const renderer = new THREE.WebGLRenderer();
renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);

// Fetch exoplanet data
fetch('https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=SELECT * FROM ps$exoplanets')
  .then(response => response.json())
  .then(data => {
    data.forEach(exoplanet => {
      const geometry = new THREE.SphereGeometry(0.1, 32, 32);
      const material = new THREE.MeshBasicMaterial({ color: 0x0077ff });
      const exoplanetMesh = new THREE.Mesh(geometry, material);
      
      // Example positions (you'll need to calculate based on real data)
      exoplanetMesh.position.set(exoplanet.distance, 0, 0);
      scene.add(exoplanetMesh);
    });
  });

// Animation loop
function animate() {
  requestAnimationFrame(animate);
  renderer.render(scene, camera);
}
animate();
